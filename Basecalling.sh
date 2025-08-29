#!/usr/bin/env bash
#SBATCH --partition=gpu_p
#SBATCH --qos=gpu_normal
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --job-name=dorado_demux_trim_chopper


# --- user config ---
DORADO_BIN="/home/haicu/dorado-0.7.1-linux-x64/bin/dorado"
MODEL="/home/haicu/ela.sauerborn/dorado-0.7.1-linux-x64/bin/dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
KIT="SQK-RBK114-24"
INPUT_DIR="./run/"   # directory containing POD5/FAST5
OUT="./"                                   # output root

# Conda envs
CONDA_ASSEMBLY_ENV="assembly"       # samtools
CONDA_PREPROC_ENV="preprocessing"   # chopper

# --- dirs ---
DEMUX_DIR="$OUT/demux"
SUMMARY_DIR="$OUT/summary"
QC_DIR="$OUT/qc_q10_l1000"
mkdir -p "$DEMUX_DIR" "$SUMMARY_DIR" "$QC_DIR"

RUN_ID="$(basename "${INPUT_DIR%/}")"
BAM="$OUT/${RUN_ID}.bam"

echo "[INFO] Basecalling -> BAM (no trim so demux sees barcodes)"
"$DORADO_BIN" basecaller \
  --device cuda:all \
  --emit-moves \
  --no-trim \
  "$MODEL" \
  --recursive "$INPUT_DIR" \
  > "$BAM" 2> "$SUMMARY_DIR/basecaller.log"

[[ -s "$BAM" ]] || { echo "[ERROR] Dorado did not produce BAM: $BAM"; exit 1; }

echo "[INFO] Demultiplex BAM -> per-barcode BAMs"
"$DORADO_BIN" demux \
  --kit-name "$KIT" \
  --output-dir "$DEMUX_DIR" \
  "$BAM" \
  > "$SUMMARY_DIR/demux.log" 2>&1

echo "[INFO] Trim per-barcode BAM -> FASTQ -> Chopper (Q>=10, len>=1000)"
# DNA-only trimming: dorado trim | samtools fastq | chopper
for bbam in "$DEMUX_DIR"/barcode*/reads.bam "$DEMUX_DIR"/barcode*.bam; do
  [[ -f "$bbam" ]] || continue
  bdir="$(basename "$(dirname "$bbam")")"
  bname="$(basename "${bbam%.*}")"
  barcode="$([[ $bdir == barcode* ]] && echo "$bdir" || echo "$bname")"
  outfq="$QC_DIR/${barcode}.q10l1000.fastq"

  "$DORADO_BIN" trim "$bbam" \
  | conda run -n "$CONDA_ASSEMBLY_ENV" samtools fastq - \
  | conda run -n "$CONDA_PREPROC_ENV" chopper -q 10 -l 1000 \
  > "$outfq"

  echo "[DONE] $outfq"
done

echo "[ALL DONE] $(date)"


