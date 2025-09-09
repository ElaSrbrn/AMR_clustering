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
DORADO_BIN="/lustre/groups/hpc/urban_lab/tools/dorado-0.9.1-linux-x64/bin/"
READS_DIR="/lustre/groups/hpc/urban_lab/backup/plasmid_project/work_package01/CSF_2025/20250226_1401_P2S-01622-B_PAY37630_aa601ade/pod5_skip/"
KIT_NAME="SQK-RBK114-96"
OUTDIR="./basecalled_v5_4mC5mC"

# Models: simplex (base) + mod model
BASE_MODEL="/lustre/groups/hpc/urban_lab/tools/dorado-0.9.1-linux-x64/bin/dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
MOD_4MC5MC="/lustre/groups/hpc/urban_lab/tools/dorado-0.9.1-linux-x64/bin/dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v3"

mkdir -p "$OUTDIR"

$DORADO_BIN/dorado basecaller sup,6mA,4mC_5mC -r $READS_DIR --models-directory /lustre/groups/hpc/urban_lab/tools/dorado-0.9.1-linux-x64/bin --kit-name $KIT_NAME > $OUTDIR/all.bam



mkdir -p "$OUTDIR"

echo "[1/3] Basecalling (4mC+5mC)…"
"$DORADO_BIN" basecaller \
  --kit-name "$KIT_NAME" \
  --trim all \
  --modified-bases-models "$MOD_4MC5MC" \
  "$BASE_MODEL" \
  -r "$READS_DIR" > "$OUTDIR/all.bam"

echo "[2/3] Demultiplexing…"
mkdir -p "$OUTDIR/demux_bam"
"$DORADO_BIN" demux \
  --kit-name "$KIT_NAME" \
  --output-dir "$OUTDIR/demux_bam" \
  "$OUTDIR/all.bam"

if ! compgen -G "$OUTDIR/demux_bam/*.bam" > /dev/null; then
  echo "ERROR: No demuxed BAMs in $OUTDIR/demux_bam (check KIT_NAME and basecalling)."
  exit 1
fi

echo "[3/3] BAM → FASTQ.GZ…"
eval "$(micromamba shell hook --shell=bash)"; micromamba activate assembly
mkdir -p "$OUTDIR/demux_fastq"
for b in "$OUTDIR"/demux_bam/*.bam; do
  base=$(basename "${b%.bam}")
  samtools fastq "$b" | gzip > "$OUTDIR/demux_fastq/${base}.fastq.gz"
done

echo "✅ Done: $OUTDIR"


