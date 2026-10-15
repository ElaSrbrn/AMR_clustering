#!/bin/bash
#SBATCH --job-name=pcs_reads
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --exclude=cpusrv08
#SBATCH --output=logs/12_read_support_%j.out
# Map 60x reads to the reference plasmid and the Flye 60x contig; report regions with <25% of median depth.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source config.sh
make_jobtmp
eval "$(micromamba shell hook --shell bash)"
micromamba activate "$ENV_MAP"
OUT=work/read_support; mkdir -p $OUT
for PID in "$@"; do
  ID=${PID%%_*}
  READS=$(ls $BASE/reads_60x_S1/*isolate_${ID}_*.fastq* 2>/dev/null | head -1)
  [ -z "$READS" ] && { echo "$PID: no 60x reads found"; continue; }
  REF=work/compare/$PID/ref.fasta
  FLYE=$(awk -F'\t' -v p="$PID" '$1==p && $6=="flye" && $7=="60x" {print $9}' results/stability_long.tsv)
  FA=$(awk -F'\t' -v i="$ID" '$1==i && $2=="flye" && $3=="60x" {print $4}' work/manifest.tsv)
  : > $JOBTMP/flye_contig.fasta
  [ -n "$FLYE" ] && [ -n "$FA" ] && $PY - "$FA" "$FLYE" > $JOBTMP/flye_contig.fasta <<'PY'
import sys; sys.path.insert(0, "lib")
from common import read_fasta
for n, h, s in read_fasta(sys.argv[1]):
    if n == sys.argv[2]:
        print(">flye_" + n); print(s)
PY
  cat $REF $JOBTMP/flye_contig.fasta > $OUT/${PID}_targets.fasta
  minimap2 -ax map-ont -t 8 --secondary=no $OUT/${PID}_targets.fasta "$READS" 2>/dev/null \
    | samtools sort -@ 4 -o $OUT/${PID}.bam - && samtools index $OUT/${PID}.bam
  samtools depth -a $OUT/${PID}.bam > $OUT/${PID}.depth
  echo "== $PID  (reads: $(basename "$READS"))"
  $PY - $OUT/${PID}.depth <<'PY'
import sys, statistics
from collections import defaultdict
dep = defaultdict(list)
for l in open(sys.argv[1]):
    c, p, d = l.split(); dep[c].append(int(d))
for c, v in dep.items():
    med = statistics.median(v); thr = 0.25 * med
    low, start = [], None
    for i, d in enumerate(v + [med]):
        if d < thr and start is None: start = i
        elif d >= thr and start is not None:
            if i - start >= 500: low.append((start + 1, i, round(statistics.mean(v[start:i]), 1)))
            start = None
    print(f"  {c}: length {len(v)}  median depth {med}  low-depth regions (<25% of median, >=500 bp): "
          + ("; ".join(f"{a}-{b} (mean {m})" for a, b, m in low) or "none"))
PY
done
