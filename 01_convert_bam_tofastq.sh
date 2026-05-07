#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# Config
THREADS=${THREADS:-8}
INDIR=${INDIR:-./}                 # BAMs: *_barcodeXX.bam
FASTQ_DIR=${FASTQ_DIR:-../fastq}
MIN_BARCODE=${MIN_BARCODE:-1}
MAX_BARCODE=${MAX_BARCODE:-24}

mkdir -p "$FASTQ_DIR"

# Convert BAM → FASTQ.gz
for bam in "$INDIR"/*_barcode*.bam; do
  barcode=$(basename "$bam" .bam | grep -oE 'barcode[0-9]{2}')
  [[ -z "$barcode" ]] && continue
  num=${barcode#barcode}

  # numeric check (handle leading zeros correctly)
  if (( 10#$num < MIN_BARCODE || 10#$num > MAX_BARCODE )); then
    continue
  fi

  base="$(basename "$bam" .bam)"
  out="$FASTQ_DIR/${base}.fastq.gz"
  [[ -s "$out" ]] && continue

  echo "Converting $bam → $out"
  samtools fastq -@ "$THREADS" -F 0x900 "$bam" | gzip > "$out"
done

#only use primary mapped reads
