#!/usr/bin/env bash
# bart_array.sh  --  run with bash on the login node, NOT with sbatch.
# Submits one BART task (MLST + AMR) per FASTQ in the given read directories.
# Outputs: <reads_dir>/bart/<sample>/{<sample>_mlst.tsv,<sample>_amr.tsv,run_info.txt}
#
# Usage:
#   bash bart_array.sh [-m scheme_map.tsv] [-p max_parallel] <reads_dir> [<reads_dir> ...]
#   e.g. bash bart_array.sh -p 20 carblot_pos/10x carblot_pos/20x carblot_neg/10x carblot_neg/20x
# scheme_map.tsv (optional, tab-separated): <sample prefix> <BART scheme name from 'bart-update -s'>

set -euo pipefail

schememap=""; max_parallel=20
while getopts "m:p:" o; do
    case $o in m) schememap="$(cd "$(dirname "$OPTARG")" && pwd)/$(basename "$OPTARG")";; p) max_parallel=$OPTARG;; *) exit 1;; esac
done
shift $((OPTIND-1))
(( $# > 0 )) || { echo "Usage: bash $0 [-m scheme_map.tsv] [-p max_parallel] <reads_dir> ..." >&2; exit 1; }
[[ -z "${SLURM_JOB_ID:-}" ]] || { echo "ERROR: run this with bash, not sbatch" >&2; exit 1; }

job="$(cd "$(dirname "$0")" && pwd)/bart_job.sh"
[[ -f "$job" ]] || { echo "ERROR: missing $job" >&2; exit 1; }
BART_ENV="${BART_ENV:-bart}"

# --- Databases must exist before tasks start (compute nodes may lack internet;
#     parallel tasks would otherwise all try to build them at once) -------------
db=$(micromamba run -n "$BART_ENV" python -c 'import bart,os;print(os.path.dirname(bart.__file__)+"/db")')
n_idx=$(ls "$db/indexes" 2>/dev/null | grep -vc '^amr' || true)
if (( n_idx == 0 )); then
    echo "ERROR: no MLST indexes in $db. Build once on the login node:  micromamba run -n $BART_ENV bart-update -p" >&2; exit 1
fi
if ! ls "$db/indexes" | grep -q '^amr'; then
    echo "ERROR: no AMR index in $db. Build once on the login node:  micromamba run -n $BART_ENV bart-update -amr" >&2; exit 1
fi

stamp="$(date +%Y%m%d_%H%M%S)"
logdir="$PWD/logs_bart_${stamp}"
mkdir -p "$logdir"
list="$logdir/reads.list"; : > "$list"

shopt -s nullglob
for d in "$@"; do
    [[ -d "$d" ]] || { echo "ERROR: not a directory: $d" >&2; exit 1; }
    d="$(cd "$d" && pwd)"
    group="$(basename "$(dirname "$d")")/$(basename "$d")"      # e.g. carblot_pos/10x
    for fq in "$d"/*.fastq "$d"/*.fastq.gz "$d"/*.fq "$d"/*.fq.gz; do
        s="$(basename "$fq")"; s="${s%.gz}"; s="${s%.fastq}"; s="${s%.fq}"
        printf '%s\t%s\t%s\t%s\n' "$group" "$s" "$fq" "$d/bart/$s" >> "$list"
    done
done
shopt -u nullglob

n=$(wc -l < "$list")
(( n > 0 )) || { echo "ERROR: no FASTQ files found" >&2; exit 1; }
dups=$(cut -f1,2 "$list" | sort | uniq -d)
[[ -z "$dups" ]] || { echo "ERROR: duplicate sample names within a group:" >&2; echo "$dups" >&2; exit 1; }
cut -f1 "$list" | sort | uniq -c
awk -F'\t' '{print NR "\t" $1 "\t" $2}' "$list" > "$logdir/task_map.tsv"

{
    echo "date: $(date -Iseconds)"
    echo "bart_db: $db"
    echo "mlst_indexes: $n_idx"
    echo "amr_index_files: $(ls -l --time-style=+%F "$db"/indexes/amr* | awk '{print $6, $7}' | tr '\n' ' ')"
} > "$logdir/bart_db_info.txt"

jobid=$(sbatch --parsable --array="1-${n}%${max_parallel}" \
    -o "$logdir/%A_%a.out" -e "$logdir/%A_%a.err" "$job" "$list" "$schememap")
echo "$jobid" > "$logdir/jobid.txt"
echo "Submitted array job $jobid: $n samples, max $max_parallel running"
echo "Not completed: sacct -j $jobid -X --format=JobID,State,ExitCode,Elapsed,MaxRSS | grep -v COMPLETED"
echo "Combine:       python3 bart_combine.py bart_results $*"
