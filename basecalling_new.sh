#!/bin/bash
#SBATCH -p gpu_p
#SBATCH -q gpu_long
#SBATCH --mem=250G
#SBATCH -t 72:00:00
#SBATCH --nice=10000
#SBATCH --job-name=dorado_nanoplot
#SBATCH -c 12
# Set working directory to dorado v1.1.1 binary location
cd /home/haicu/ela.sauerborn/dorado-1.1.1-linux-x64/bin

# === Input / Output Paths ===
READS_DIR="/lustre/groups/hpc/urban_lab/backup/plasmid_project/work_package01/CSF_2025/96Pool/20250226_1401_P2S-01622-B_PAY37630_aa601ade/pod5_skip"
OUTDIR="/lustre/groups/hpc/urban_lab/projects/plasmid_project/stored_resistant_enterobacterales/work_package01/26022025/basecalled_v5_6mA4mC5mC"
DEMUX="${OUTDIR}/demultiplexed"
TRIMDIR="${OUTDIR}/trimmed"
mkdir -p "$OUTDIR" "$DEMUX" "$TRIMDIR"

# === Models ===
BASE_MODEL="dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
MODBASE_MODELS="dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v3,dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v3"
KIT="SQK-RBK114-96"

# === Step 1: Basecalling with modbase support ===
./dorado basecaller "$BASE_MODEL" "$READS_DIR" \
  --kit-name "$KIT" \
  --modified-bases-models "$MODBASE_MODELS" \
  --output "$OUTDIR/all.bam" \
  --verbose

# === Step 2: Demultiplexing ===
./dorado demux --kit-name "$KIT" \
  --output-dir "$DEMUX" \
  "$OUTDIR/all.bam"

# === Step 3: Trimming ===
for bam in "$DEMUX"/*.bam; do
    sample="$(basename "$bam" .bam)"
    mkdir -p "${TRIMDIR}/${sample}"
    ./dorado trim "$bam" --kit-name "$KIT" > "${TRIMDIR}/${sample}/trimmed.bam"
done



eval "$(micromamba shell hook --shell=bash)"; 
micromamba activate assembly

mkdir -p "$OUTDIR/demux_fastq"
for b in "$OUTDIR"/demux_bam/*.bam; do
  base=$(basename "${b%.bam}")
  samtools fastq "$b" | gzip > "$OUTDIR/demux_fastq/${base}.fastq.gz"
done



