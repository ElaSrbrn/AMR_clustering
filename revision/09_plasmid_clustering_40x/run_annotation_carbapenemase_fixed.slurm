#!/usr/bin/env bash
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --nice=10000
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL
#SBATCH --job-name=annot_carb
#SBATCH --output=annot_%x_%j.out
#SBATCH --error=annot_%x_%j.err

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch $0 <polished_fasta_dir>" >&2
  echo "Example: sbatch $0 polished_fasta_40x" >&2
  exit 1
fi

POLISHED_DIR="$1"
if [[ ! -d "$POLISHED_DIR" ]]; then
  echo "ERROR: input directory not found: $POLISHED_DIR" >&2
  exit 1
fi

THREADS="${SLURM_CPUS_PER_TASK:-16}"
base="$(basename "$POLISHED_DIR")"
coverage="$(printf '%s\n' "$base" | grep -Eo '[0-9]+x' | head -n1 || true)"
if [[ -z "$coverage" ]]; then
  coverage="$base"
fi

OUTROOT="./annotation_${coverage}"
SUMMARY_DIR="${OUTROOT}/summaries"
CARB_PLASMID_DIR="${OUTROOT}/carbapenemase_encoding_plasmids"
mkdir -p "$OUTROOT" "$SUMMARY_DIR" "$CARB_PLASMID_DIR"

eval "$(micromamba shell hook --shell=bash)"

mapfile -t ASMS < <(
  find "$POLISHED_DIR" -maxdepth 1 -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \) | sort
)

echo "Found ${#ASMS[@]} polished assemblies in ${POLISHED_DIR}"
if [[ ${#ASMS[@]} -eq 0 ]]; then
  echo "ERROR: no polished assemblies found in ${POLISHED_DIR}" >&2
  exit 1
fi

for asm in "${ASMS[@]}"; do
  [[ -e "$asm" ]] || continue

  fname="$(basename "$asm")"
  sample="${fname%.fasta}"
  sample="${sample%.fa}"
  sample="${sample%.fna}"

  outdir="${OUTROOT}/${sample}"
  amrdir="${outdir}/amrfinder"
  mobdir="${outdir}/mob_recon"

  mkdir -p "$outdir" "$amrdir" "$mobdir"

  echo "==> $sample"

  micromamba run -n amr_env amrfinder \
    -n "$asm" \
    --threads "$THREADS" \
    > "${amrdir}/${sample}_assembly_amrfinder.tsv" \
    2> "${amrdir}/${sample}_assembly_amrfinder.log" || \
    echo "AMRFinder failed on assembly for $sample" >> "${outdir}/pipeline_warnings.txt"

  micromamba run -n mobsuite_env mob_recon \
    --infile "$asm" \
    --outdir "$mobdir" \
    --force \
    2>> "${outdir}/pipeline_warnings.txt" || \
    echo "MOB-recon failed for $sample" >> "${outdir}/pipeline_warnings.txt"

done

python3 /mnt/data/build_carbapenemase_summary.py \
  --annotation-root "$OUTROOT" \
  --assembly-root "$POLISHED_DIR" \
  --coverage "$coverage"

echo
echo "Done."
echo "Annotation root:                 $OUTROOT"
echo "Summary dir:                     $SUMMARY_DIR"
echo "Carbapenemase plasmid dir:       $CARB_PLASMID_DIR"
