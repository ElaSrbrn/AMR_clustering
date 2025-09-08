#!/bin/bash

# root with polished assemblies and *_mobsuite outputs
ROOT="assemblies/polished"
OUT="assemblies/plasmids_by_name"
mkdir -p "$OUT"

# optional index of what was collected
INDEX="$OUT/index.tsv"
echo -e "plasmid\tbarcode\tsource_path" > "$INDEX"

# collect plasmid*.fa*
while IFS= read -r -d '' f; do
  mobsuite_dir="$(dirname "$f")"                          # .../<sample>_mobsuite
  sample_dir="$(basename "$mobsuite_dir" "_mobsuite")"    # strip suffix → e.g. barcode0102
  base="$(basename "$f")"

  # accept plasmidXX.fa / plasmidXX.fasta
  if [[ "$base" =~ ^plasmid([A-Za-z0-9_.-]+)\.fa(sta)?$ ]]; then
    pid="${BASH_REMATCH[1]}"
  else
    echo "[WARN] skipping non-plasmid file: $f" >&2
    continue
  fi

  outdir="$OUT/plasmid${pid}"
  mkdir -p "$outdir"
  outfa="$outdir/${sample_dir}.fasta"

  cp -f "$f" "$outfa"
  echo -e "plasmid${pid}\t${sample_dir}\t${f}" >> "$INDEX"
done < <(find "$ROOT" -type f -path "*/plasmid*.fa*" -print0 | sort -zV)

echo "[OK] grouped plasmids → $OUT"
echo "Index written: $INDEX"
