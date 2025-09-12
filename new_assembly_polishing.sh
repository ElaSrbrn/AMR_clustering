#!/bin/bash
#SBATCH -p gpu_p
#SBATCH -q gpu
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH --qos gpu_normal
#SBATCH --partition=gpu_p
#SBATCH --nice=10000
#SBATCH --mail-user=ela.sauerborn@helmholtz-munich.de
#SBATCH --mail-type=ALL
#SBATCH --gres=gpu:2
#SBATCH --job-name=assembly_new
#SBATCH -c 12

MERGED_FASTQ="/lustre/groups/hpc/urban_lab/projects/plasmid_project/stored_resistant_enterobacterales/work_package01/02052025/csf_repeat_02052025/basecalled_v5_6mA4mC5mC/filtered_fastq"
# Where to store unpolished Flye outputs
ASM_UNPOLISHED="/lustre/groups/hpc/urban_lab/projects/plasmid_project/stored_resistant_enterobacterales/work_package01/02052025/csf_repeat_02052025/basecalled_v5_6mA4mC5mC/assemblies/unpolished"

# Where to store polished Medaka outputs
ASM_POLISHED="/lustre/groups/hpc/urban_lab/projects/plasmid_project/stored_resistant_enterobacterales/work_package01/02052025/csf_repeat_02052025/basecalled_v5_6mA4mC5mC/assemblies/polished"

# Make sure they exist
mkdir -p "$ASM_UNPOLISHED" "$ASM_POLISHED"

# ========================
# Step 1: Assembly with Flye
# ========================
eval "$(micromamba shell hook --shell=bash)"
micromamba activate assembly

for fq in "$MERGED_FASTQ"/barcode*_filtered.fastq; do
  [[ -e "$fq" ]] || continue
  sample=$(basename "${fq%_filtered.fastq}")

  draft_dir="$ASM_UNPOLISHED/${sample}_flye"
  draft="$draft_dir/assembly.fasta"
  polished_dir="$ASM_POLISHED/${sample}_medaka"
  polished="$polished_dir/consensus.fasta"

  if [[ ! -s "$draft" ]]; then
    mkdir -p "$draft_dir"
    flye \
      --nano-hq "$fq" \
      --out-dir "$draft_dir" \
      --genome-size 5m \
      --asm-coverage 100
  fi
done

micromamba deactivate

# ========================
# Step 2: Polishing with Medaka
# ========================
micromamba activate medaka_env

for fq in "$MERGED_FASTQ"/barcode*_filtered.fastq; do
  [[ -e "$fq" ]] || continue
  sample=$(basename "${fq%_filtered.fastq}")

  draft_dir="$ASM_UNPOLISHED/${sample}_flye"
  draft="$draft_dir/assembly.fasta"
  polished_dir="$ASM_POLISHED/${sample}_medaka"
  polished="$polished_dir/consensus.fasta"

  if [[ ! -s "$polished" ]]; then
    minimap2 -x map-ont -t "$THREADS" "$draft" "$fq" \
      | samtools sort -@ "$THREADS" -o "${draft%.fasta}.reads.bam"
    samtools index "${draft%.fasta}.reads.bam"

    mkdir -p "$polished_dir"
    medaka consensus \
      -i "$fq" \
      -d "$draft" \
      -o "$polished_dir" \
      -t 12 \
      --bacteria

    rm -f "${draft%.fasta}.reads.bam" "${draft%.fasta}.reads.bam.bai"
  fi

  echo "→ ${sample}: unpolished=$draft polished=$polished"
done

echo "All done."
