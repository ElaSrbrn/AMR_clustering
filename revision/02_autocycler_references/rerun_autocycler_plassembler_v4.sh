#!/bin/bash
#SBATCH --job-name=ac_plass_v4
#SBATCH --exclude=cpusrv08
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs_ac_rerun4/%x_%A_%a.out
#
# Adds Plassembler assemblies to existing Autocycler runs and rebuilds the consensus.
# The subsampled reads were deleted after the original run, so they are regenerated first.
# autocycler subsample uses seed 0 by default, so the same input + settings should give identical
# subsamples; this is verified against the original subsample.yaml and the isolate is skipped if not.
#
# Usage (from the revision_scripts folder):
#   1) Check that every isolate's 60x read file is found (no SLURM needed, quick):
#        MODE=check bash rerun_autocycler_plassembler_v2.sh
#   2) Test one isolate:
#        mkdir -p logs_ac_rerun4
#        sbatch --array=1 rerun_autocycler_plassembler_v2.sh
#   3) Run the rest:
#        sbatch --array=2-$(wc -l < isolate_dirs.txt)%20 rerun_autocycler_plassembler_v2.sh
#   4) Summary:
#        grep -h "^STATUS" logs_ac_rerun4/*.out | cut -f3 | sort | uniq -c
#        grep -h "^STATUS" logs_ac_rerun4/*.out | grep -v OK

set -uo pipefail

# Per-job temporary directory on node-local disk (the default /tmp filled up on some nodes)
if [ -d /localscratch ] && [ -w /localscratch ]; then TMPBASE=/localscratch/$USER; else TMPBASE=$(pwd)/tmp_jobs; fi
JOBTMP=$TMPBASE/${SLURM_JOB_ID:-x}_${SLURM_ARRAY_TASK_ID:-0}
mkdir -p "$JOBTMP"
export TMPDIR=$JOBTMP TMP=$JOBTMP TEMP=$JOBTMP
trap 'rm -rf "$JOBTMP"' EXIT
echo "TMPDIR: $JOBTMP ($(df -h "$JOBTMP" | awk 'NR==2{print $4}') free)"
LIST=${LIST:-isolate_dirs.txt}

# Locate the 60x read file for an isolate folder .../60x/autocycler_fast/<name>
# by looking in .../60x/ for <name>.fastq(.gz) or files starting with <name> or <name without _60x>.
find_reads() {
    local dir=$1 name base readsdir
    name=$(basename "$dir")
    base=${name%_60x}
    readsdir=$(dirname "$(dirname "$dir")")
    local c=() f
    for f in "$readsdir/$name".fastq "$readsdir/$name".fastq.gz "$readsdir/$name".fq "$readsdir/$name".fq.gz; do
        [ -f "$f" ] && c+=( "$f" )
    done
    shopt -s nullglob
    if [ ${#c[@]} -eq 0 ]; then c=( "$readsdir/$name"*.f*q* ); fi
    if [ ${#c[@]} -eq 0 ]; then c=( "$readsdir/$base"*.f*q* ); fi
    shopt -u nullglob
    if [ ${#c[@]} -eq 1 ]; then echo "${c[0]}"; return 0; fi
    echo "${#c[@]} candidates in $readsdir: ${c[*]:-none}"
    return 1
}

if [ "${MODE:-run}" = "check" ]; then
    while read -r d; do
        if r=$(find_reads "$d"); then
            echo -e "OK\t$(basename "$d")\t$r"
        else
            echo -e "PROBLEM\t$(basename "$d")\t$r"
        fi
    done < "$LIST"
    exit 0
fi

eval "$(micromamba shell hook --shell bash)"
AC_ENV=${AC_ENV:-autocycler_0.8.0}
micromamba activate "$AC_ENV"
export PLASSEMBLER_DB=${PLASSEMBLER_DB:-/ictstr01/home/haicu/ela.sauerborn/mamba/envs/autocycler/plassembler_db}
echo "Autocycler: $(autocycler --version)  env: $AC_ENV  PLASSEMBLER_DB: $PLASSEMBLER_DB"

DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
NAME=$(basename "$DIR")
status() { echo -e "STATUS\t$NAME\t$1"; }
echo "Isolate folder: $DIR"
cd "$DIR" || { status FAIL_cd; exit 1; }

READS=$(find_reads "$DIR") || { echo "$READS"; status FAIL_reads_not_found; exit 1; }
echo "Reads: $READS"

JOBLOG=assemblies/joblog.tsv
GS=$(grep -o -- '--genome_size [0-9]*' "$JOBLOG" | head -1 | awk '{print $2}')
[ -n "$GS" ] || { status FAIL_no_genome_size; exit 1; }
ORIG_YAML=subsampled_reads/subsample.yaml
[ -f "$ORIG_YAML" ] || { status FAIL_no_original_yaml; exit 1; }

# 1. Regenerate the subsamples into a separate folder and verify they match the originals
rm -rf subsampled_reads_regen
autocycler subsample --reads "$READS" --out_dir subsampled_reads_regen --genome_size "$GS" --count 2 \
    || { status FAIL_subsample; exit 1; }
if ! diff -q "$ORIG_YAML" subsampled_reads_regen/subsample.yaml > /dev/null; then
    echo "Original vs regenerated subsample.yaml:"
    diff "$ORIG_YAML" subsampled_reads_regen/subsample.yaml
    rm -f subsampled_reads_regen/*.fastq
    status FAIL_subsample_mismatch
    exit 1
fi
echo "Subsamples identical to the original run."

# 2. Run Plassembler with exactly the original options, pointed at the regenerated subsamples
rm -f assemblies/plassembler_*.fasta
PL_FAIL=0
while read -r cmd; do
    cmd=${cmd//subsampled_reads\//subsampled_reads_regen\/}
    echo "Running: $cmd"
    eval "$cmd" || { echo "Plassembler returned non-zero"; PL_FAIL=1; }
done < <(grep -o 'autocycler helper plassembler.*' "$JOBLOG")

rm -f subsampled_reads_regen/*.fastq   # keep only the yaml, as in the original run

if [ "$PL_FAIL" -eq 1 ]; then
    rm -f assemblies/plassembler_*.fasta
    status FAIL_plassembler   # consensus left untouched
    exit 1
fi
find assemblies -maxdepth 1 -name 'plassembler_*.fasta' -empty -delete
N_PL=$(ls assemblies/plassembler_*.fasta 2>/dev/null | wc -l)
echo "Non-empty Plassembler assemblies: $N_PL"
if [ "$N_PL" -eq 0 ]; then
    status OK_no_small_plasmids_consensus_unchanged
    exit 0
fi

# 3. Rebuild the consensus. The original consensus (no Plassembler) is kept in autocycler_out_noplassembler.
#    If that folder already exists, autocycler_out is the redundant rebuild from the first attempt and is replaced.
if [ -d autocycler_out_noplassembler ]; then
    rm -rf autocycler_out
elif [ -d autocycler_out ]; then
    mv autocycler_out autocycler_out_noplassembler
fi

autocycler compress -i assemblies -a autocycler_out            || { status FAIL_compress; exit 1; }
autocycler cluster -a autocycler_out                           || { status FAIL_cluster; exit 1; }
for c in autocycler_out/clustering/qc_pass/cluster_*; do
    autocycler trim -c "$c"    || { status FAIL_trim; exit 1; }
    autocycler resolve -c "$c" || { status FAIL_resolve; exit 1; }
done
autocycler combine -a autocycler_out -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa \
    || { status FAIL_combine; exit 1; }

status "OK_rebuilt_with_${N_PL}_plassembler"
