#!/usr/bin/env bash
# autocycler_fast_array.sh  --  run with bash on the login node, NOT with sbatch.
# Writes a list of all FASTQ files in a directory and submits autocycler_fast_job.sh
# as a SLURM array (one isolate per task). Resources come from the #SBATCH header
# of autocycler_fast_job.sh.
#
# Usage:
#   bash autocycler_fast_array.sh <fastq_dir> [max_parallel]
# Run it from the directory that contains autocycler_fast_job.sh and autocycler_fast.sh
# (outputs go to ./autocycler_fast/<isolate>/).

set -euo pipefail

fastq_dir="${1:?Usage: bash $0 <fastq_dir> [max_parallel]}"
max_parallel="${2:-4}"
job="$PWD/autocycler_fast_job.sh"

[[ -d "$fastq_dir" ]]           || { echo "ERROR: not a directory: $fastq_dir" >&2; exit 1; }
[[ -f "$job" ]]                 || { echo "ERROR: missing $job" >&2; exit 1; }
[[ -f "$PWD/autocycler_fast.sh" ]] || { echo "ERROR: missing $PWD/autocycler_fast.sh" >&2; exit 1; }
[[ -z "${SLURM_JOB_ID:-}" ]]    || { echo "ERROR: run this with bash, not sbatch" >&2; exit 1; }

fastq_dir="$(cd "$fastq_dir" && pwd)"
stamp="$(date +%Y%m%d_%H%M%S)"
logdir="$PWD/logs_autocycler_${stamp}"
mkdir -p "$logdir"
list="$logdir/reads.list"

shopt -s nullglob
files=( "$fastq_dir"/*.fastq "$fastq_dir"/*.fastq.gz )
shopt -u nullglob
(( ${#files[@]} > 0 )) || { echo "ERROR: no .fastq/.fastq.gz files in $fastq_dir" >&2; exit 1; }
printf '%s\n' "${files[@]}" | sort > "$list"
n=$(wc -l < "$list")

# Output folders are named after the file stem, so stems must be unique
dups=$(sed -E 's#.*/##; s/\.fastq(\.gz)?$//' "$list" | sort | uniq -d)
[[ -z "$dups" ]] || { echo "ERROR: same isolate name twice (fastq and fastq.gz?):" >&2; echo "$dups" >&2; exit 1; }

awk '{print NR "\t" $0}' "$list" > "$logdir/task_map.tsv"

jobid=$(sbatch --parsable \
    --array="1-${n}%${max_parallel}" \
    -o "$logdir/%A_%a.out" -e "$logdir/%A_%a.err" \
    "$job" "$list")

echo "$jobid" > "$logdir/jobid.txt"
echo "Submitted array job $jobid: $n isolates, max $max_parallel running at once"
echo "Task -> FASTQ map: $logdir/task_map.tsv"
echo "Status:            squeue -j $jobid"
echo "Not completed:     sacct -j $jobid -X --format=JobID,State,ExitCode,Elapsed,MaxRSS | grep -v COMPLETED"
echo "Re-run task 7:     sbatch --array=7 -o $logdir/%A_%a.out -e $logdir/%A_%a.err $job $list"
