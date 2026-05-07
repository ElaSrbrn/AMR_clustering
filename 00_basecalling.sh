#!/bin/bash
#SBATCH -p gpu_p
#SBATCH -q gpu_normal
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH --mail-user=samir.vargas@helmholtz-munich.de
#SBATCH --mail-type=ALL
#SBATCH --gres=gpu:1
#SBATCH --job-name=dorado
#SBATCH --cpus-per-task=8


READS_DIR="pod5"
OUTDIR="basecalled_v5allmods"
DEMUX="${OUTDIR}/demultiplexed"
DORADO_BIN="/home/haicu/ela.sauerborn/dorado-1.1.1-linux-x64/bin/dorado"
mkdir -p "$DEMUX" #"$TRIMDIR"
# BASE_MODEL="dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
# ODBASE_MODELS="dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v3,dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v3"
KIT="SQK-RBK114-96"

#basecalling
${DORADO_BIN} basecaller sup@v5.0.0,4mC_5mC,6mA $READS_DIR --kit-name $KIT --no-trim > "$OUTDIR/all.bam"

#demultiplexing
${DORADO_BIN} demux --kit-name "$KIT" --output-dir "$DEMUX" "$OUTDIR/all.bam"
