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
#SBATCH --job-name=dorado_nanoplot
#SBATCH -c 2

CONFIG_FILE="/lustre/groups/hpc/urban_lab/tools/dorado-0.9.1-linux-x64/bin/dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v3"
READS_DIR="/lustre/groups/hpc/urban_lab/backup/plasmid_project/work_package01/CSF_2025/20250226_1401_P2S-01622-B_PAY37630_aa601ade"   # directory with raw ONT reads
KIT_NAME="SQK-RBK114-96"                 # your kit
DORADO_BIN="/lustre/groups/hpc/urban_lab/tools/dorado-0.9.1-linux-x64/bin/dorado"
OUTDIR="./basecalled_v5_mod"

mkdir -p "$OUTDIR"

"$DORADO_BIN" basecaller \
  --kit-name "$KIT_NAME" \
  --trim \
  --emit-bam \
  "$CONFIG_FILE" \
  -r "$READS_DIR" > "$OUTDIR/all.bam"

mkdir -p "$OUTDIR/demux_bam"
"$DORADO_BIN" demux \
  --kit-name "$KIT_NAME" \
  --output-dir "$OUTDIR/demux_bam" \
  "$OUTDIR/all.bam"

# Activate micromamba for samtools
eval "$(micromamba shell hook --shell=bash)"
micromamba activate assembly

mkdir -p "$OUTDIR/demux_fastq"
for b in "$OUTDIR"/demux_bam/*.bam; do
  base=$(basename "${b%.bam}")
  samtools fastq "$b" | gzip > "$OUTDIR/demux_fastq/${base}.fastq.gz"
done
