#!/usr/bin/env bash
# autocycler_fast_submit.sh
# Builds a list of FASTQ files and submits autocycler_fast.sh as a SLURM array
# (one task per FASTQ). autocycler_fast.sh must read line ${SLURM_ARRAY_TASK_ID}
# of the list file passed as $1 (as in its existing input block).
#
# Usage:
#   bash autocycler_fast_submit.sh <fastq_dir> [max_parallel] [script]
# Example:
#   bash autocycler_fast_submit.sh /lustre/.../reads 10 ./autocycler_fast.sh

set -euo pipefail
#activate autocycler_0.8.0
# --- Arguments ---------------------------------------------------------------
fastq_dir="${1:?Usage: $0 <fastq_dir> [max_parallel] [script]}"
max_parallel="${2:-10}"
script="${3:-./autocycler_fast.sh}"

# --- Resources per task (edit as needed) -------------------------------------
PARTITION="cpu_p"
QOS="cpu_normal"
CPUS=16
MEM="64G"
TIME="12:00:00"
NICE=10000
MAIL=""

# --- Checks ------------------------------------------------------------------
[[ -d "$fastq_dir" ]] || { echo "ERROR: not a directory: $fastq_dir" >&2; exit 1; }
[[ -f "$script" ]]    || { echo "ERROR: script not found: $script" >&2; exit 1; }
command -v sbatch >/dev/null || { echo "ERROR: sbatch not found (run on a login node)" >&2; exit 1; }

fastq_dir="$(cd "$fastq_dir" && pwd)"     # absolute path, so tasks find the files
script="$(cd "$(dirname "$script")" && pwd)/$(basename "$script")"

# --- Build the list ----------------------------------------------------------
stamp="$(date +%Y%m%d_%H%M%S)"
logdir="$(pwd)/logs_autocycler_${stamp}"
mkdir -p "$logdir"
list="${logdir}/reads.list"

shopt -s nullglob
files=( "$fastq_dir"/*.fastq "$fastq_dir"/*.fastq.gz "$fastq_dir"/*.fq "$fastq_dir"/*.fq.gz )
shopt -u nullglob
(( ${#files[@]} > 0 )) || { echo "ERROR: no FASTQ files in $fastq_dir" >&2; exit 1; }

printf '%s\n' "${files[@]}" | sort > "$list"
n=$(wc -l < "$list")

# Same sample present twice (e.g. x.fastq and x.fastq.gz)?
dups=$(sed -E 's#.*/##; s/\.(fastq|fq)(\.gz)?$//' "$list" | sort | uniq -d)
if [[ -n "$dups" ]]; then
    echo "ERROR: duplicate sample names (compressed and uncompressed copies?):" >&2
    echo "$dups" >&2
    exit 1
fi

# Task -> sample map, for the record and for re-running single tasks
awk '{print NR "\t" $0}' "$list" > "${logdir}/task_map.tsv"

echo "FASTQ dir:   $fastq_dir"
echo "Files:       $n"
echo "Max running: $max_parallel"
echo "Script:      $script"
echo "List:        $list"
echo "Logs:        $logdir"
(( n <= 98 )) || echo "NOTE: $n files found (more than 98); all will be submitted." >&2

# --- Submit ------------------------------------------------------------------
jobid=$(sbatch --parsable \
    -J autocycler \
    -p "$PARTITION" -q "$QOS" \
    -c "$CPUS" --mem="$MEM" -t "$TIME" \
    --nice="$NICE" \
    --mail-type=END,FAIL --mail-user="$MAIL" \
    --array="1-${n}%${max_parallel}" \
    -o "${logdir}/%A_%a.out" -e "${logdir}/%A_%a.err" \
    "$script" "$list")

echo "Submitted array job ${jobid} (tasks 1-${n})"
echo "$jobid" > "${logdir}/jobid.txt"
echo
echo "Check status:      squeue -j ${jobid}"
echo "Failed tasks:      sacct -j ${jobid} --format=JobID,State,ExitCode,Elapsed,MaxRSS | grep -v COMPLETED"
echo "Re-run task 7:     sbatch --array=7 [same options] $script $list"
