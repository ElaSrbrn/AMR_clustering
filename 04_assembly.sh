#!/usr/bin/env bash
#SBATCH --mem=100G
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=16
#SBATCH --nice=10000

#to run this script, you have to specify the coverages as follows: sbatch 04_assembly.sh 10x (until 100x) manually
# there is also an option to parallelize the assembly to speed things up by splitting the assembly into arrays, meaning 
#that fastq files are assembled in parallel. 


COV="${1:?Usage: sbatch scripts/01_flye.sh <coverage_dir>}"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
GENOME_SIZE="${GENOME_SIZE:-5m}"

READS_DIR="./${COV}"
ASM_UNPOLISHED="./flye_${COV}"

mkdir -p "$ASM_UNPOLISHED"

[[ -d "$READS_DIR" ]] || { echo "ERROR: missing $READS_DIR" >&2; exit 1; }

eval "$(micromamba shell hook --shell=bash)"
micromamba activate assembly

mapfile -t FASTQS < <(
  find "$READS_DIR" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) | sort
)

[[ ${#FASTQS[@]} -gt 0 ]] || { echo "ERROR: no FASTQs in $READS_DIR" >&2; exit 1; }

for fq in "${FASTQS[@]}"; do
  sample="$(basename "$fq")"
  sample="${sample%.fastq.gz}"
  sample="${sample%.fastq}"

  draft_dir="${ASM_UNPOLISHED}/${sample}_flye"
  draft="${draft_dir}/assembly.fasta"

  if [[ -s "$draft" ]]; then
    echo "SKIP Flye: $sample"
    continue
  fi

  echo "RUN Flye: $sample"
  flye \
    --nano-hq "$fq" \
    --out-dir "$draft_dir" \
    --genome-size "$GENOME_SIZE" \
    --threads "$THREADS"
done
