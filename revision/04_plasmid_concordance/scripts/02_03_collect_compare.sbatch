#!/bin/bash
#SBATCH --job-name=pcs_compare
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --exclude=cpusrv08
#SBATCH --output=logs/02_03_compare_%j.out
#
# Steps 02 + 03 + 04: reference plasmids, pairing, alignment/classification, summaries.
# Usage (from the pipeline folder, after all 01 tasks report n/n annotated):
#   sbatch scripts/02_03_collect_compare.sbatch
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source config.sh
make_jobtmp
export LEN_TOL MIN_AF MIN_ID MIN_BLOCK BREAK_GAP CHROM_MIN CHIMERA_UNALIGNED FRAGMENT_MIN_BP SMALL_PLASMID
$PY scripts/02_collect.py   || exit 1
$PY scripts/03_compare.py   || exit 1
$PY scripts/04_summarise.py || exit 1
