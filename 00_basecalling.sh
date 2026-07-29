#!/usr/bin/env bash
set -euo pipefail
#SBATCH -p gpu_p
#SBATCH -q gpu_normal
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --gres=gpu:1
#SBATCH --job-name=dorado
#SBATCH --cpus-per-task=8


READS_DIR="pod5"
OUTDIR="basecalled_v5allmods"
DEMUX="${OUTDIR}/demultiplexed"
# Path to the dorado binary. Uses whatever is on PATH by default; override with
#   DORADO_BIN=/path/to/dorado sbatch 00_basecalling.sh
DORADO_BIN="${DORADO_BIN:-dorado}"
mkdir -p "$DEMUX" #"$TRIMDIR"
# BASE_MODEL="dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
# ODBASE_MODELS="dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v3,dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v3"
KIT="SQK-RBK114-96"

# Basecalling with the super-accuracy model. The 4mC_5mC and 6mA modified-base
# models are also called, so the output BAM carries methylation tags. These are
# not used downstream in this pipeline; drop the two suffixes from the model
# string below if they are not wanted.
${DORADO_BIN} basecaller sup@v5.0.0,4mC_5mC,6mA $READS_DIR --kit-name $KIT --no-trim > "$OUTDIR/all.bam"

#demultiplexing
${DORADO_BIN} demux --kit-name "$KIT" --output-dir "$DEMUX" "$OUTDIR/all.bam"
