#!/usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH --nice=10000
#SBATCH --job-name=bart

# Read-level MLST and AMR gene screening with BART (one FASTQ per array task).
# Do not sbatch directly; use:  bash bart_array.sh <reads_dir> [<reads_dir> ...]
# $1 = list file: <group>\t<sample>\t<fastq>\t<outdir>     $2 = optional scheme map

set -eo pipefail
BART_ENV="${BART_ENV:-bart}"
AMR_ID="${AMR_ID:-90}"            # BART -amr percent identity (BART default 90)
eval "$(micromamba shell hook --shell=bash)"
micromamba activate "$BART_ENV"
set -u

list="${1:?list file required}"
schememap="${2:-}"
[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo "ERROR: run via bart_array.sh" >&2; exit 1; }

IFS=$'\t' read -r group sample fq outdir < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$list")
[[ -f "$fq" ]] || { echo "ERROR: reads not found: $fq" >&2; exit 1; }
echo "Task ${SLURM_ARRAY_TASK_ID}: $group / $sample"

THREADS="${SLURM_CPUS_PER_TASK:-8}"
mkdir -p "$outdir"
# Private temp dir: BART writes <sample>.* and choose_scheme.bart into $TMPDIR
export TMPDIR="$outdir/tmp_bart"
rm -rf "$TMPDIR"; mkdir -p "$TMPDIR"

# Optional forced scheme: longest prefix in the map that matches the sample name
scheme=""
if [[ -n "$schememap" && -f "$schememap" ]]; then
    scheme=$(awk -F'\t' -v s="$sample" '
        $1!="" && (s==$1 || index(s, $1"_")==1) && length($1)>best {best=length($1); sc=$2}
        END {print sc}' "$schememap")
fi
sarg=(); [[ -n "$scheme" ]] && sarg=(-s "$scheme")
echo "MLST scheme: ${scheme:-chosen by BART (refseq_masher)}"

# -k on both runs: without it BART crashes in its cleanup step (non-f-string path bug)
# A failure in one run must not prevent the other: record exit codes instead of aborting.
mlst_rc=0; amr_rc=0
bart "$fq" -r ont -t "$THREADS" -vv -k "${sarg[@]}" \
    > "$outdir/${sample}_mlst.tsv" 2> "$outdir/${sample}_mlst.log" || mlst_rc=$?

bart "$fq" -r ont -t "$THREADS" -amr "$AMR_ID" -k \
    > "$outdir/${sample}_amr.tsv" 2> "$outdir/${sample}_amr.log" || amr_rc=$?
echo "exit codes: mlst=$mlst_rc amr=$amr_rc"
(( mlst_rc == 0 )) || { echo "MLST run failed, last log lines:" >&2; tail -5 "$outdir/${sample}_mlst.log" >&2; }
(( amr_rc == 0 ))  || { echo "AMR run failed, last log lines:" >&2;  tail -5 "$outdir/${sample}_amr.log" >&2; }

{
    echo "group: $group"
    echo "sample: $sample"
    echo "reads: $fq"
    echo "bart: $(grep -m1 -o 'this is bart [^ ]*' "$outdir/${sample}_mlst.log" || echo unknown)"
    echo "kma: $(kma -v 2>&1 | head -1)"
    echo "mlst_scheme_forced: ${scheme:-no}"
    echo "amr_percid: $AMR_ID"
    echo "exit_codes: mlst=$mlst_rc amr=$amr_rc"
    echo "date: $(date -Iseconds)"
} > "$outdir/run_info.txt"

rm -rf "$TMPDIR"
echo "Done: $outdir"
(( mlst_rc == 0 && amr_rc == 0 ))   # task shows FAILED in sacct if either run failed
