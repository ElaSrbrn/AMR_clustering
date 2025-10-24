#!/bin/bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH --nice=10000
#SBATCH --mail-user=ela.sauerborn@helmholtz-munich.de
#SBATCH --mail-type=ALL
#SBATCH --job-name=chopper

INDIR="./fastq"
OUTDIR="./trimmed_fastq"

mkdir -p "$OUTDIR"

for f in "$INDIR"/*barcode*.fastq.gz; do
  [ -e "$f" ] || continue
  base=$(basename "$f" .fastq.gz)

  echo "Processing: $base"

  # Filter with Chopper (q=9, l=1000) and remove duplicate read IDs (keep first)
  chopper -q 9 -l 1000 -i "$f" | \
    awk 'NR%4==1{if(seen[$1]++)skip=1;else skip=0}!skip' \
    > "$OUTDIR/${base}_filtered_q9l1000.fastq"
done
