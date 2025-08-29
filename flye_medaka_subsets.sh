#!/usr/bin/env bash
#SBATCH --partition=cpu_p
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --job-name=flye_medaka_subsets_bacteria
#SBATCH --mail-type=ALL

# -------- user config --------
SEQTK_BIN="/home/haicu/ela.sauerborn/seqtk/seqtk"   
INPUT_FASTQ="/Chopper_output_path" #update
GENOME_SIZE=5000000
DEPTHS="10 20 30 40 50 60 70 80 90 100"
SEED=42
ENV_NAME="assembly"                                 
OUT="./assemble_subsets"

# Activate env 
eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"

# Dirs
READS_TAG="$(basename "${INPUT_FASTQ%.*}")"
SUBDIR="$OUT/${READS_TAG}"
SHUF_DIR="$SUBDIR/shuffle"
SUBSETS_DIR="$SUBDIR/subsets"
FLYE_DIR="$SUBDIR/flye"
MED_DIR="$SUBDIR/medaka"

mkdir -p "$SHUF_DIR" "$SUBSETS_DIR" "$FLYE_DIR" "$MED_DIR"


SHUF_FQ="$SHUF_DIR/${READS_TAG}.shuffled.fastq"
"$SEQTK_BIN" sample -s "$SEED" "$INPUT_FASTQ" 1 > "$SHUF_FQ"

TARGETS_TSV="$SHUF_DIR/targets.tsv"; : > "$TARGETS_TSV"
for d in $DEPTHS; do echo -e "${d}\t$((GENOME_SIZE * d))" >> "$TARGETS_TSV"; done

COUNTS_TSV="$SHUF_DIR/depth_read_counts.tsv"
awk -v OFS='\t' '
  NR==FNR { tgt[$1]=$2; depths[n++]=$1; next }
  (NR%4)==2 { cum+=length($0); reads+=1; while(i<n && cum>=tgt[depths[i]]){print depths[i],reads; i++} }
  END { if(i<n){ printf("[WARN] Insufficient total bases for:"); for(;i<n;i++) printf(" %sx",depths[i]); printf("\n") > "/dev/stderr" } }
' "$TARGETS_TSV" "$SHUF_FQ" > "$COUNTS_TSV"

while IFS=$'\t' read -r depth reads_needed; do
  lines=$((reads_needed * 4))
  subset="$SUBSETS_DIR/${READS_TAG}.depth_${depth}x.fastq"
  echo "  - ${depth}x → ${reads_needed} reads"
  head -n "$lines" "$SHUF_FQ" > "$subset"
done < "$COUNTS_TSV"
ln -sf "$SHUF_FQ" "$SUBSETS_DIR/${READS_TAG}.depth_full.fastq"

for subset in "$SUBSETS_DIR"/*.fastq; do
  base="$(basename "$subset" .fastq)"
  outdir="$FLYE_DIR/$base"; mkdir -p "$outdir"
  flye --nano-hq "$subset" \
       --genome-size "$GENOME_SIZE" \
       --out-dir "$outdir" \

done

for subset in "$SUBSETS_DIR"/*.fastq; do
  base="$(basename "$subset" .fastq)"
  draft="$FLYE_DIR/$base/assembly.fasta"
  mdir="$MED_DIR/$base"; mkdir -p "$mdir"
  medaka_consensus \
    -i "$subset" \
    -d "$draft" \
    -o "$mdir" \
    --bacteria \
done


