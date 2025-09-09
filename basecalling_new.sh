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
#SBATCH --job-name=dorado_nanoplot
#SBATCH -c 2

# Paths (adjust only these)
DORADO_BIN="path_to_urban_lab/tools/dorado-0.9.1-linux-x64/bin/"
READS_DIR="/your_path/pod5_skip/"
KIT_NAME="SQK-RBK114-96"
OUTDIR="./basecalled_v5_mod"

mkdir -p "$OUTDIR"

$DORADO_BIN/dorado basecaller sup,6mA,4mC_5mC -r $READS_DIR --models-directory /path_to_urban_lab/tools/dorado-0.9.1-linux-x64/bin --kit-name $KIT_NAME > $OUTDIR/all.bam


mkdir -p "$OUTDIR/demux_bam"
"$DORADO_BIN" demux \
  --kit-name "$KIT_NAME" \
  --output-dir "$OUTDIR/demux_bam" \
  "$OUTDIR/all.bam"


eval "$(micromamba shell hook --shell=bash)"; 
micromamba activate assembly

mkdir -p "$OUTDIR/demux_fastq"
for b in "$OUTDIR"/demux_bam/*.bam; do
  base=$(basename "${b%.bam}")
  samtools fastq "$b" | gzip > "$OUTDIR/demux_fastq/${base}.fastq.gz"
done



