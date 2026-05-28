#!/usr/bin/env bash

# Config
THREADS=${THREADS:-8}
INDIR=${INDIR:-./}                 
FASTQ_DIR=${FASTQ_DIR:-../fastq}
MIN_BARCODE=${MIN_BARCODE:-1}
MAX_BARCODE=${MAX_BARCODE:-24}

for bam in "$INDIR"/*_barcode*.bam; do
  base=${bam##*/}
  base=${base%.bam}

  [[ $base =~ barcode([0-9]{2}) ]] || continue
  num=$((10#${BASH_REMATCH[1]}))

  (( num < MIN_BARCODE || num > MAX_BARCODE )) && continue

  out="$FASTQ_DIR/${base}.fastq.gz"
  [[ -s $out ]] || samtools fastq -@ "$THREADS" -F 0x900 "$bam" | gzip > "$out"
done
