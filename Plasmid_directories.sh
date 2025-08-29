#!/usr/bin/env bash
#SBATCH --partition=cpu_p
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --job-name=cluster_plasmids_by_name


# Root produced by your earlier pipeline
ROOT="./assemble_subsets"
MOB_DIR="$ROOT/reports/mobsuite"
OUT="$ROOT/reports/plasmid_name_clusters"
mkdir -p "$OUT"

# 1) Collect plasmid FASTAs (filenames like plasmidA455.fasta, plasmid_1.fa, etc.)
mapfile -t PFASTAS < <(find "$MOB_DIR" -type f \
  \( -iname 'plasmid*.fa' -o -iname 'plasmid*.fasta' \) \
  | sort)

if (( ${#PFASTAS[@]} == 0 )); then
  echo "[ERROR] No plasmid FASTAs under $MOB_DIR"; exit 1
fi

# 2) Build an index: name \t subset \t path
INDEX="$OUT/plasmids_index_by_name.tsv"
: > "$INDEX"
for fa in "${PFASTAS[@]}"; do
  name="$(basename "$fa")"                       # e.g., plasmidA455.fasta
  subset="$(basename "$(dirname "$fa")")"        # e.g., 0102_deduped.depth_40x
  printf "%s\t%s\t%s\n" "$name" "$subset" "$fa" >> "$INDEX"
done

# 3) Names seen >=2 times
NAMES_MULTI="$OUT/names_multi.txt"
awk -F'\t' '{c[$1]++} END{for(n in c) if(c[n]>=2) print n}' "$INDEX" | sort > "$NAMES_MULTI"

if ! [[ -s "$NAMES_MULTI" ]]; then
  echo "[INFO] No plasmid filenames shared across outputs."; exit 0
fi

# 4) For each shared name, make a dir and symlink members as <subset>__<name>
while IFS= read -r name; do
  cdir="$OUT/$name"
  mkdir -p "$cdir"

  # manifest
  awk -F'\t' -v N="$name" 'BEGIN{OFS="\t"; print "name","subset","path"} $1==N{print $1,$2,$3}' \
    "$INDEX" > "$cdir/cluster.tsv"

  # symlinks
  while IFS=$'\t' read -r nm subset path; do
    target="$(readlink -f "$path" 2>/dev/null || realpath "$path" 2>/dev/null || echo "$path")"
    ln -sf "$target" "$cdir/${subset}__${nm}"
  done < <(awk -F'\t' -v N="$name" '$1==N{print $1"\t"$2"\t"$3}' "$INDEX")
done < "$NAMES_MULTI"

# 5) Global summary
awk -F'\t' 'NR==FNR{hit[$1]=1; next} ($1 in hit){print}' "$NAMES_MULTI" "$INDEX" \
| sort -k1,1 -k2,2 \
> "$OUT/shared_plasmid_names.tsv"

echo "[DONE] Clusters in: $OUT"
echo "  - Per-name dirs: $OUT/<plasmidName>/{cluster.tsv, <subset>__<plasmidName>}"
echo "  - Summary     : $OUT/shared_plasmid_names.tsv"
