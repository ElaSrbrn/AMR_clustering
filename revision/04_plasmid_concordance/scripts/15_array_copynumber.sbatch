#!/bin/bash
#SBATCH --job-name=pcs_array
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --exclude=cpusrv08
#SBATCH --output=logs/15_array_copynumber_%j.out
# Step 15: copy number of the 9,335 bp IS26-IS6100 tandem unit of 7365_2 per spanning read.
# Usage: sbatch scripts/15_array_copynumber.sbatch   or   READS=/path/to/full_reads.fastq.gz sbatch scripts/15_array_copynumber.sbatch
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source ./config.sh
make_jobtmp
eval "$(micromamba shell hook --shell bash)"
micromamba activate "$ENV_MAP"
PID=${PID:-7365_2}; ID=${PID%%_*}
LA=${LA:-26000}; RA=${RA:-64000}; UNIT=${UNIT:-9335}; REF_UNITS=${REF_UNITS:-3}
READS=${READS:-$(ls $BASE/reads_60x_S1/*isolate_${ID}_*.fastq* 2>/dev/null | head -1)}
[ -s "$READS" ] || { echo "no reads for $ID"; exit 1; }
OUT=work/array_copynumber/$PID; mkdir -p $OUT
echo "== $PID  reads: $READS"
CAT=cat; case "$READS" in *.gz) CAT="zcat";; esac
echo "reads >= 40 kb in the read set: $($CAT "$READS" | awk 'NR%4==2 && length($0)>=40000{n++} END{print n+0}')"
minimap2 -c -x map-ont -t 8 --secondary=no work/compare/$PID/ref.fasta "$READS" 2>/dev/null > $OUT/reads.paf
$PY scripts/15_array_copynumber.py $OUT/reads.paf $LA $RA $UNIT $REF_UNITS $OUT/per_read.tsv
cp $OUT/per_read.tsv results/array_copynumber_${PID}.tsv
echo "written: results/array_copynumber_${PID}.tsv"
