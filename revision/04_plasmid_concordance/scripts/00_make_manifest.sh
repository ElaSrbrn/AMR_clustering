#!/bin/bash
# Step 00: build work/manifest.tsv with one row per isolate x assembler x depth and the assembly path.
# Usage (from the pipeline folder): bash scripts/00_make_manifest.sh
source "$(dirname "$0")/../config.sh"
cd "$PCS" || exit 1
mkdir -p work results logs
OUT=work/manifest.tsv
printf 'id\tassembler\tdepth\tpath\tstatus\n' > $OUT
shopt -s nullglob

override () {  # id assembler depth -> path from overrides.tsv, if any
  awk -F'\t' -v i="$1" -v a="$2" -v d="$3" 'NR>1 && $1==i && $2==a && $3==d {print $4}' overrides.tsv
}
emit () {  # id assembler depth candidates...
  local id=$1 asm=$2 dep=$3; shift 3
  local ov; ov=$(override "$id" "$asm" "$dep")
  if [ -n "$ov" ]; then set -- "$ov"; fi
  if   [ $# -eq 0 ]; then printf '%s\t%s\t%s\t\tMISSING\n' "$id" "$asm" "$dep" >> $OUT
  elif [ $# -gt 1 ]; then printf '%s\t%s\t%s\t%s\tMULTIPLE\n' "$id" "$asm" "$dep" "$*" >> $OUT
  elif [ ! -s "$1" ]; then printf '%s\t%s\t%s\t%s\tEMPTY\n' "$id" "$asm" "$dep" "$1" >> $OUT
  else printf '%s\t%s\t%s\t%s\tOK\n' "$id" "$asm" "$dep" "$1" >> $OUT
  fi
}

while read -r ID; do
  [ -z "$ID" ] && continue
  # Reference: Autocycler 60x consensus
  dir=$(grep "/isolate_${ID}_" "$ISOLATE_DIRS" | head -1)
  c=(); [ -n "$dir" ] && [ -s "$dir/autocycler_out/consensus_assembly.fasta" ] && c=("$dir/autocycler_out/consensus_assembly.fasta")
  emit "$ID" autocycler 60x "${c[@]}"
  for D in $DEPTHS; do
    # Hybracter (complete or incomplete)
    emit "$ID" hybracter "$D" $HYB_OUT/${ID}_${D}/FINAL_OUTPUT/*/${ID}_${D}_final.fasta
    # Flye
    g=${FLYE_GLOB//\{ID\}/$ID}; g=${g//\{DEPTH\}/$D}
    emit "$ID" flye "$D" $g
  done
done < "$IDS_FILE"

echo "Manifest: $PCS/$OUT"
awk -F'\t' 'NR>1{print $2, $3, $5}' $OUT | sort | uniq -c
echo "Rows needing attention (MISSING / MULTIPLE / EMPTY):"
awk -F'\t' 'NR>1 && $5!="OK"' $OUT | head -40
