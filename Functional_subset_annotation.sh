#!/usr/bin/env bash
#SBATCH --partition=cpu_p
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --job-name=amr_mob_on_subsets
#SBATCH --mail-type=ALL
set -euo pipefail

ROOT="./assemble_subsets"
AMR_DIR="$ROOT/reports/amrfinder"
MOB_DIR="$ROOT/reports/mobsuite"
mkdir -p "$AMR_DIR" "$MOB_DIR"

# Find all Medaka-polished assemblies
mapfile -t FASTAS < <(find "$ROOT" -type f \
  \( -name 'consensus.fasta' -o -name 'consensus.fa' -o -name 'consensus_polished.fasta' \) \
  -path '*/medaka/*/*' | sort)

if (( ${#FASTAS[@]} == 0 )); then
  echo "[ERROR] No polished consensus FASTA found under $ROOT (expected at */medaka/<subset>/consensus*.fasta)"
  exit 1
fi

THREADS="${SLURM_CPUS_PER_TASK:-1}"

for fa in "${FASTAS[@]}"; do
  subset="$(basename "$(dirname "$fa")")"
  echo "[INFO] Subset: $subset"

  # --- AMRFinderPlus ---
  amr_out="$AMR_DIR/${subset}.amrfinder.tsv"
  amr_log="$AMR_DIR/${subset}.log"
  conda run -n amrfinder amrfinder \
    -n "$fa" \
    --plus \
    --organism Bacteria \
    --threads "$THREADS" \
    -o "$amr_out" \
    >"$amr_log" 2>&1
  echo "  AMR -> $amr_out"

  # --- MOB-suite ---
  mob_out="$MOB_DIR/${subset}"
  mkdir -p "$mob_out"
  mob_log="$mob_out/mob_recon.log"
  conda run -n mobsuite mob_recon \
    -i "$fa" \
    -o "$mob_out" \
    -n "$THREADS" \
    >"$mob_log" 2>&1
  echo "  MOB -> $mob_out"
done


