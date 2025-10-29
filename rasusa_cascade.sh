#!/usr/bin/env bash
set -euo pipefail

# ----- EDIT THESE -----
G="5m"          # Genome size (e.g., 5m)
SEED=42         # Random seed
OUTDIR_BASE="../"  # Parent directory for output folders
# ----------------------

mkdir -p "${OUTDIR_BASE}/100x" "${OUTDIR_BASE}/70x" "${OUTDIR_BASE}/40x" "${OUTDIR_BASE}/10x"

normalize_name() {
  # $1 = basename; $2 = cov tag (100x/70x/40x/10x)
  local bn="$1" cov="$2"
  local stem="${bn%.fastq}"
  if [[ "$stem" == *"_filtered" ]]; then
    stem="${stem%_filtered}"
  elif [[ "$stem" == *"_merged" ]]; then
    stem="${stem%_merged}"
  fi
  printf "%s_%s.fastq" "$stem" "$cov"
}

shopt -s nullglob
for f in *.fastq; do
  bn="$(basename "$f")"

  out100="$(normalize_name "$bn" 100x)"
  out70="$(normalize_name "$bn" 70x)"
  out40="$(normalize_name "$bn" 40x)"
  out10="$(normalize_name "$bn" 10x)"

  echo "Processing $bn ..."

  rasusa reads -g "$G" -c 100 -s "$SEED" "$f" > "${OUTDIR_BASE}/100x/$out100"
  rasusa reads -g "$G" -c 70  -s "$SEED" "${OUTDIR_BASE}/100x/$out100" > "${OUTDIR_BASE}/70x/$out70"
  rasusa reads -g "$G" -c 40  -s "$SEED" "${OUTDIR_BASE}/70x/$out70"  > "${OUTDIR_BASE}/40x/$out40"
  rasusa reads -g "$G" -c 10  -s "$SEED" "${OUTDIR_BASE}/40x/$out40"  > "${OUTDIR_BASE}/10x/$out10"
done

