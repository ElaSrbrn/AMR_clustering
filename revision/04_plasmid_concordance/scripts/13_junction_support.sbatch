#!/bin/bash
#SBATCH --job-name=pcs_junct
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --exclude=cpusrv08
#SBATCH --output=logs/13_junction_support_%j.out
# Step 13: junction-spanning reads (see header comments in the attached file for details)
# Usage: sbatch scripts/13_junction_support.sbatch 6057_2 7419_3 7365_2 7475_3
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source config.sh
make_jobtmp
eval "$(micromamba shell hook --shell bash)"
micromamba activate "$ENV_MAP"
export W=${W:-10000} OVERHANG=${OVERHANG:-3000} MINID=${MINID:-0.80} BIGINDEL=${BIGINDEL:-100}
OUT=work/junction_support; mkdir -p $OUT results
PIDS=("$@"); [ ${#PIDS[@]} -gt 0 ] || PIDS=(6057_2 7419_3 7365_2 7475_3)
: > $OUT/all_rows.tsv
for PID in "${PIDS[@]}"; do
  ID=${PID%%_*}; D=$OUT/$PID; rm -rf $D; mkdir -p $D/win $D/paf
  READS=$(ls $BASE/reads_60x_S1/*isolate_${ID}_*.fastq* 2>/dev/null | head -1)
  [ -z "$READS" ] && { echo "$PID: no 60x reads"; continue; }
  echo "== $PID  (reads: $(basename "$READS"))"
  $PY scripts/13_junction_windows.py "$PID" "$D" || continue
  for FA in $D/win/*.fasta; do
    K=$(basename $FA .fasta)
    minimap2 -c -x map-ont -t 8 --secondary=no "$FA" "$READS" 2>/dev/null > $D/paf/$K.paf
  done
  $PY scripts/13_junction_count.py "$D" >> $OUT/all_rows.tsv
done
{ printf "plasmid_id\tflye_depth\tbreakpoint\tclass\tkind\tside\tref_pos\tquery_pos\treads_ref_only\treads_flye_only\treads_both\tfraction_flye\tinterpretation\n"; cat $OUT/all_rows.tsv; } > results/junction_support.tsv
echo "written: results/junction_support.tsv  (W=$W, overhang=$OVERHANG, min identity=$MINID, indel cut-off=$BIGINDEL)"
column -t -s $'\t' results/junction_support.tsv
