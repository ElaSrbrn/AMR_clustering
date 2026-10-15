#!/bin/bash
#SBATCH --job-name=hybracter_medaka
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=logs_hybracter_medaka/%x_%A_%a.out
#
# One Hybracter long-single run per row of hybracter_depth_series.csv (isolate x depth),
# with Hybracter's default Medaka polishing, using the Hybracter container (database and all sub-environments are inside the image,
# so nothing is downloaded on the offline compute nodes and no -d is needed).
#
# Usage (from hybracter_work):
#   mkdir -p logs_hybracter_medaka
#   sbatch --array=1 run_hybracter_array.sh                                   # test one sample
#   sbatch --array=2-$(wc -l < hybracter_depth_series.csv)%20 run_hybracter_array.sh

set -uo pipefail
module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true
RUN=$(command -v apptainer || command -v singularity) || { echo "no apptainer/singularity found"; exit 1; }
export PYTHONNOUSERSITE=1

SIF=/ictstr01/groups/hpc/urban_lab/projects/plasmid_project/containers/hybracter.sif
SHEET=${SHEET:-/ictstr01/groups/hpc/urban_lab/projects/plasmid_project/revision_scripts/hybracter_depth_series.csv}
OUT=/ictstr01/groups/hpc/urban_lab/projects/plasmid_project/hybracter_depth_out_medaka
mkdir -p "$OUT"

IFS=, read -r SAMPLE READS MINLEN < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SHEET")
echo "Sample: $SAMPLE  Reads: $READS  Min chromosome: $MINLEN"
READS=$(readlink -f "$READS")      # resolve the symlink so the real file is visible inside the container
[ -s "$READS" ] || { echo -e "STATUS\t$SAMPLE\tFAIL_reads"; exit 1; }

mkdir -p "$OUT/$SAMPLE"
cd "$OUT/$SAMPLE"      # one working dir per sample, so parallel Snakemake runs do not share .snakemake/
"$RUN" exec --bind /ictstr01 "$SIF" \
  hybracter long-single -l "$READS" -s "$SAMPLE" -c "$MINLEN" -o "$OUT/$SAMPLE" -t 16 \
  --skip_qc --circular_chromosome
rc=$?
if [ $rc -eq 0 ]; then echo -e "STATUS\t$SAMPLE\tOK"; else echo -e "STATUS\t$SAMPLE\tFAIL_hybracter_rc$rc"; fi
exit $rc
