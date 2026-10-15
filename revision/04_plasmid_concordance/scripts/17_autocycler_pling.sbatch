#!/bin/bash
#SBATCH --job-name=pcs_ac_pling
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/17_autocycler_pling_%j.out
# Step 17b: Mash (<=0.005) and pling (containment 0.3, DCJ-Indel <=4) on the 60x Autocycler reference plasmids,
# same script and settings as the Flye clustering; then pairwise comparison with Flye 60x.
# Usage: PLING_ENV=<env name> sbatch scripts/17_autocycler_pling.sbatch
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source ./config.sh
: "${PLING_ENV:?set PLING_ENV to the environment that provides mash and pling}"
NEWCRE=/ictstr01/groups/hpc/urban_lab/projects/plasmid_project/stored_resistant_enterobacterales/new_cre
GROUPPY=$NEWCRE/carblot_pos/scripts/group_plasmids_mash_then_pling_v2_compat.py
FLYE60=$NEWCRE/carblot_pos/annotation_60x
eval "$(micromamba shell hook --shell bash)"
micromamba activate "$PLING_ENV"
echo "mash $(mash --version)"; pling --version 2>/dev/null || pling --help 2>&1 | head -2
$PY scripts/17_autocycler_clustering_input.py
cd clustering_autocycler60x
python3 "$GROUPPY" --annotation-glob "annotation_autocycler_60x" --mash-threshold 0.005 --pling-containment 0.3 --threads 8
cat annotation_autocycler_60x/mash_then_pling/group_summary.tsv 2>/dev/null | cut -f1-8
ln -sfn "$FLYE60" annotation_flye_60x
python3 ../scripts/11_compare_pling_depths.py annotation_flye_60x annotation_autocycler_60x IncN IncL_M
