#!/bin/bash
#SBATCH --job-name=pcs_context
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/16_context_synteny_%j.out
# Step 16: synteny of the carbapenemase genetic context (gene +/- 10 kb) between reference and test assemblies.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source ./config.sh
export TMPDIR=$PWD/work/tmp16; mkdir -p $TMPDIR
export CONDA_PREFIX=${MAMBA_ROOT_PREFIX:-$HOME/mamba}/envs/amr_env
ROOT=${MAMBA_ROOT_PREFIX:-$HOME/mamba}
export PATH=$ROOT/envs/$ENV_MAP/bin:$ROOT/envs/amr_env/bin:$PATH
echo "minimap2 $(minimap2 --version)"; amrfinder --version
$PY scripts/16_context_synteny.py "$@"
