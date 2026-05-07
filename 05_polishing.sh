#!/bin/bash
#SBATCH -p cpu_p
#SBATCH -q cpu
#SBATCH --mem=128G
#SBATCH -t 24:00:00


#similar to the assembly you have to specify the depth you want to polish, so sbatch 05_polishing.sh 10x (up to 100x or whatever depth you are interested in)
#Again, this could be parallelized as well to speed things up. You can also link the assembly with the polishing script to be directly added as a next step for each assembled isolate

set -euo pipefail


THREADS="${SLURM_CPUS_PER_TASK:-16}"
MEDAKA_MODEL="${MEDAKA_MODEL:-r1041_e82_400bps_bacterial_methylation}"
COPY_FINAL_FASTA="${COPY_FINAL_FASTA:-0}"

READS_DIR="./${COV}"
ASM_UNPOLISHED="./flye_${COV}"
ASM_POLISHED="./polished_flye_${COV}"
POLISHED_FASTA_DIR="./polished_fasta_${COV}"

mkdir -p "$ASM_POLISHED" "$POLISHED_FASTA_DIR"

[[ -d "$READS_DIR" ]] || { echo "ERROR: missing $READS_DIR" >&2; exit 1; }
[[ -d "$ASM_UNPOLISHED" ]] || { echo "ERROR: missing $ASM_UNPOLISHED" >&2; exit 1; }

eval "$(micromamba shell hook --shell=bash)"
micromamba activate assembly

command -v medaka_consensus >/dev/null || { echo "ERROR: medaka_consensus not found" >&2; exit 127; }
command -v minimap2 >/dev/null || { echo "ERROR: minimap2 not found; Medaka needs it" >&2; exit 127; }
command -v samtools >/dev/null || { echo "ERROR: samtools not found; Medaka needs it" >&2; exit 127; }

mapfile -t FASTQS < <(
  find "$READS_DIR" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) | sort
)

[[ ${#FASTQS[@]} -gt 0 ]] || { echo "ERROR: no FASTQs in $READS_DIR" >&2; exit 1; }

for fq in "${FASTQS[@]}"; do
  base="$(basename "$fq")"
  sample="${base%.fastq.gz}"
  sample="${sample%.fastq}"

  draft="${ASM_UNPOLISHED}/${sample}_flye/assembly.fasta"
  polished_dir="${ASM_POLISHED}/${sample}_medaka"
  polished_fasta="${polished_dir}/consensus.fasta"
  final_fasta="${POLISHED_FASTA_DIR}/${sample}.fasta"
  log="${polished_dir}/medaka.log"

  if [[ ! -s "$draft" ]]; then
    echo "WARNING: missing Flye draft for $sample"
    continue
  fi

  mkdir -p "$polished_dir"

  if [[ -s "$polished_fasta" ]]; then
    echo "SKIP Medaka: $sample"
  else
    echo "RUN Medaka: $sample"
    medaka_consensus \
      -i "$fq" \
      -d "$draft" \
      -o "$polished_dir" \
      -t "$THREADS" \
      -m "$MEDAKA_MODEL" \
      > "$log" 2>&1
  fi

  if [[ -s "$polished_fasta" ]]; then
    if [[ "$COPY_FINAL_FASTA" -eq 1 ]]; then
      cp -f "$polished_fasta" "$final_fasta"
    else
      ln -sfn "$(realpath "$polished_fasta")" "$final_fasta"
    fi
    echo "OK: $sample -> $final_fasta"
  else
    echo "FAIL: $sample; see $log"
  fi
done
