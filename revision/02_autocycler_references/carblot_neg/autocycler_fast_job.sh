#!/bin/bash
#SBATCH --mem=100G
#SBATCH -t 72:00:00
#SBATCH --cpus-per-task=48
#SBATCH --nice=10000
#SBATCH --mail-type=ALL
#SBATCH --job-name=autocycler_fast

# Runs Ryan Wick's autocycler_fast.sh (2 read subsets x 6 assemblers = 12 input assemblies)
# unmodified, taken from https://github.com/rrwick/Autocycler at tag v0.8.0
# (commit d184c3cddf7f1afa389bf430e03d3d4cef6a35a4,
#  pipelines/Automated_Autocycler_Bash_script_by_Ryan_Wick/autocycler_fast.sh).
# The fast script needs Autocycler >= v0.7.0 (it calls `autocycler combine -r ... -t ...`).
#
# Usage, single isolate:
#   sbatch autocycler_fast_job.sh /path/to/reads.fastq
# Usage, many isolates: do not sbatch this file directly, run the launcher
#   bash autocycler_fast_array.sh /path/to/fastq_dir [max_parallel]
# which writes the list and submits this script as an array (one FASTQ per task).

set -eo pipefail

# Environment with Autocycler >= 0.7.0 and the assemblers (override: sbatch --export=ALL,AC_ENV=name)
AC_ENV="${AC_ENV:-autocycler_0.8.0}"
eval "$(micromamba shell hook --shell=bash)"
micromamba activate "$AC_ENV"
set -u   # enabled after activation, since activation scripts can reference unset variables

WAY="${SLURM_SUBMIT_DIR:-$PWD}"
AUTO="${WAY}/autocycler_fast.sh"
[[ -f "$AUTO" ]] || { echo "ERROR: missing $AUTO" >&2; exit 1; }

# --- Version check -----------------------------------------------------------
ver=$(autocycler --version | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
if [[ "$(printf '%s\n' 0.7.0 "$ver" | sort -V | head -1)" != "0.7.0" ]]; then
    echo "ERROR: Autocycler $ver found; autocycler_fast.sh needs >= v0.7.0 (tested script is from v0.8.0)" >&2
    exit 1
fi
echo "Autocycler version: $ver"

# --- Temp dirs ---------------------------------------------------------------
mkdir -p /lustre/groups/hpc/urban_lab/projects/etcpp
export TMPDIR=/lustre/groups/hpc/urban_lab/projects/etcpp
export TEMP="$TMPDIR"
export TMP="$TMPDIR"

# --- Input: single FASTQ or one line of a list (array job) -------------------
arg="${1:-${WAY}/isolate_6737_FBD_PAU_PAY_40x.fastq}"
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    [[ -f "$arg" && "$arg" != *.fastq && "$arg" != *.fastq.gz ]] || \
        { echo "ERROR: array mode needs a list file (one FASTQ per line) as \$1, got: $arg" >&2; exit 1; }
    read=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$arg")
    [[ -n "$read" ]] || { echo "ERROR: no line ${SLURM_ARRAY_TASK_ID} in $arg" >&2; exit 1; }
    echo "Array task ${SLURM_ARRAY_TASK_ID}: $read"
else
    read="$arg"
fi
[[ -f "$read" ]] || { echo "ERROR: reads not found: $read" >&2; exit 1; }
# Absolute path, because the pipeline runs inside the output directory
read=$(realpath "$read")

readsbase=$(basename "$read")
readsname=${readsbase%.fastq.gz}
readsname=${readsname%.fastq}

# --- Threads ------------------------------------------------------------------
# Total threads = jobs x threads per job must not exceed the allocated CPUs.
CPUS=${SLURM_CPUS_PER_TASK:-16}
JOBS=3
THREADS=$(( CPUS / JOBS ))
(( THREADS >= 1 )) || { echo "ERROR: fewer CPUs ($CPUS) than parallel jobs ($JOBS)" >&2; exit 1; }
echo "Running $JOBS assembler jobs in parallel with $THREADS threads each ($CPUS CPUs allocated)"

# --- Run ---------------------------------------------------------------------
outdir="${WAY}/autocycler_fast/${readsname}"
mkdir -p "$outdir"
pushd "$outdir" >/dev/null

# autocycler_fast.sh <read_fastq> <threads> <jobs> [read_type]  (read_type defaults to ont_r10)
bash "$AUTO" "$read" "$THREADS" "$JOBS"

# Keep a record of how this assembly was made
{
    echo "isolate_reads: $read"
    echo "autocycler_version: $ver"
    echo "pipeline: autocycler_fast.sh (rrwick/Autocycler v0.8.0, d184c3c)"
    echo "pipeline_md5: $(md5sum "$AUTO" | cut -d' ' -f1)"
    echo "threads_per_job: $THREADS"
    echo "parallel_jobs: $JOBS"
    echo "date: $(date -Iseconds)"
} > run_info.txt

popd >/dev/null
echo "Done: ${outdir}"
