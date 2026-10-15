#!/bin/bash
#SBATCH --job-name=pcs_plsdb
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --exclude=cpusrv08
#SBATCH --output=logs/14_plsdb_%j.out
# Step 14 (R1.2): nearest PLSDB plasmid for every carbapenemase-encoding reference plasmid.
# Usage (from the pipeline folder): sbatch scripts/14_plsdb.sbatch
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source config.sh
make_jobtmp
PLDB=${PLASSEMBLER_DB:-/ictstr01/home/haicu/ela.sauerborn/mamba/envs/autocycler/plassembler_db}
TOOLS=${TOOLS_BIN:-${MAMBA_ROOT_PREFIX:-$HOME/mamba}/envs/autocycler/bin}
export PATH="$TOOLS:$PATH"
export TOPN=${TOPN:-5}
for t in mash minimap2 samtools; do command -v $t >/dev/null || echo "WARN: $t not found in $TOOLS"; done
command -v mash >/dev/null || { echo "mash is required; set TOOLS_BIN to an env bin folder that has mash"; exit 1; }
MSH=$(ls "$PLDB"/*plsdb*.msh 2>/dev/null | head -1)
META=$(ls "$PLDB"/*plsdb*.tsv 2>/dev/null | head -1)
FASTA=$(ls "$PLDB"/*plsdb*.fna "$PLDB"/*plsdb*.fasta "$PLDB"/*plsdb*.fa 2>/dev/null | head -1)
echo "PLSDB sketch: ${MSH:-none}"; echo "PLSDB metadata: ${META:-none}"; echo "PLSDB FASTA: ${FASTA:-none (coverage/identity skipped)}"
echo "mash $(mash --version)"; command -v minimap2 >/dev/null && echo "minimap2 $(minimap2 --version)"
[ -n "$MSH" ] || { echo "no PLSDB Mash sketch in $PLDB"; ls -l "$PLDB"; exit 1; }
OUT=work/plsdb; mkdir -p $OUT results
cut -f1 results/stability_long.tsv | tail -n +2 | sort -u > $OUT/plasmids.txt
: > $OUT/queries.fasta
while read -r PID; do
  f=work/compare/$PID/ref.fasta
  [ -s "$f" ] || { echo "missing $f"; continue; }
  awk -v p="$PID" 'NR==1{print ">"p; next} {print}' "$f" >> $OUT/queries.fasta
done < $OUT/plasmids.txt
echo "query plasmids: $(grep -c '>' $OUT/queries.fasta)"
mash dist -p 8 -i "$MSH" $OUT/queries.fasta > $OUT/mash_all.tsv
sort -t$'\t' -k2,2 -k3,3g $OUT/mash_all.tsv | awk -F'\t' -v n=$TOPN '{c[$2]++; if (c[$2]<=n) print}' > $OUT/mash_top.tsv
echo "top hits: $(wc -l < $OUT/mash_top.tsv)"
: > $OUT/aln.paf
if [ -n "$FASTA" ] && command -v minimap2 >/dev/null && command -v samtools >/dev/null; then
  [ -s "$FASTA.fai" ] || samtools faidx "$FASTA" 2>/dev/null || { cp "$FASTA" $JOBTMP/plsdb.fa; FASTA=$JOBTMP/plsdb.fa; samtools faidx "$FASTA"; }
  cut -f1 $OUT/mash_top.tsv | sort -u > $OUT/hit_ids.txt
  samtools faidx "$FASTA" -r $OUT/hit_ids.txt > $OUT/hits.fasta 2> $OUT/faidx.err
  minimap2 -c -x asm10 -t 8 $OUT/hits.fasta $OUT/queries.fasta 2>/dev/null > $OUT/aln.paf
fi
$PY scripts/14_plsdb_summary.py "$OUT" "${META:-}"
