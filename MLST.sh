#!/usr/bin/env bash
#SBATCH --partition=cpu_p
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --job-name=mlst_on_mobsuite_chromosomes
set -euo pipefail

# Root from your previous pipeline
ROOT="./assemble_subsets"
MOB_DIR="$ROOT/reports/mobsuite"
MLST_DIR="$ROOT/reports/mlst"
mkdir -p "$MLST_DIR"

# Find every MOB-suite chromosome FASTA
mapfile -t FASTAS < <(find "$MOB_DIR" -type f -name 'chromosome.fasta' | sort)
if (( ${#FASTAS[@]} == 0 )); then
  echo "[ERROR] No chromosome.fasta under $MOB_DIR"; exit 1
fi

# One combined TSV for all subsets
COMBINED="$MLST_DIR/mlst_all.tsv"
: > "$COMBINED"

for fa in "${FASTAS[@]}"; do
  subset="$(basename "$(dirname "$fa")")"           # e.g., reads.depth_10x
  out="$MLST_DIR/${subset}.mlst.tsv"
  echo "[INFO] MLST: $subset"
  # Prefix subset name, keep mlst’s native columns
  conda run -n mlst mlst "$fa" \
  | awk -v s="$subset" 'BEGIN{OFS="\t"} {print s,$0}' \
  | tee "$out" >> "$COMBINED"
done

echo "[DONE] Per-subset: $MLST_DIR/*.mlst.tsv"
echo "[DONE] Combined :  $COMBINED"
