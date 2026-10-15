#!/bin/bash
#SBATCH --job-name=pcs_annotate
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --exclude=cpusrv08
#SBATCH --output=logs/01_annotate_%A_%a.out
#
# Step 01: annotate every assembly of one isolate (reference + Flye + Hybracter at all depths).
#   contigs.tsv   contig, length, circular
#   amr.tsv       AMRFinderPlus (nucleotide mode)
#   mob.tsv       MOB-typer --multi on all contigs shorter than CHROM_MIN
# One array task per isolate (line of IDS_FILE). Finished assemblies are skipped, so the job can be resubmitted.
# Usage (from the pipeline folder):
#   sbatch --array=1 scripts/01_annotate.sbatch                                  # test one isolate
#   sbatch --array=2-$(wc -l < ids_s1a_60x.txt)%15 scripts/01_annotate.sbatch
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source config.sh
make_jobtmp
eval "$(micromamba shell hook --shell bash)"

ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$IDS_FILE")
echo "Isolate $ID  TMPDIR $JOBTMP"
awk -F'\t' -v i="$ID" 'NR>1 && $1==i && $5=="OK" {print $2"\t"$3"\t"$4}' work/manifest.tsv |
while IFS=$'\t' read -r ASM DEP FA; do
  OUT=work/annot/$ID/${ASM}_${DEP}
  [ -f "$OUT/DONE" ] && { echo "skip $ASM $DEP (done)"; continue; }
  mkdir -p "$OUT"
  echo "== $ASM $DEP  $FA"
  $PY -c "import sys; sys.path.insert(0,'lib'); from common import *; r=contig_table(sys.argv[1]); write_tsv(sys.argv[2], r, ['contig','length','circular'])" "$FA" "$OUT/contigs.tsv" \
    || { echo "FAIL contigs $ASM $DEP"; continue; }

  micromamba activate "$ENV_AMR"
  amrfinder -n "$FA" --threads 8 -o "$OUT/amr.tsv" > "$OUT/amr.log" 2>&1 || { echo "FAIL amrfinder $ASM $DEP"; continue; }

  micromamba activate "$ENV_MOB"
  awk -F'\t' -v m="$CHROM_MIN" 'NR>1 && $2<m {print $1}' "$OUT/contigs.tsv" > "$OUT/nonchrom.ids"
  if [ -s "$OUT/nonchrom.ids" ]; then
    awk 'NR==FNR{k[$1]=1; next} /^>/{split(substr($0,2),a," "); p=(a[1] in k)} p' "$OUT/nonchrom.ids" "$FA" > "$JOBTMP/nonchrom.fasta"
    mob_typer --multi -i "$JOBTMP/nonchrom.fasta" -o "$OUT/mob.tsv" -n 8 > "$OUT/mob.log" 2>&1 \
      || { echo "FAIL mob_typer $ASM $DEP"; continue; }
  else
    : > "$OUT/mob.tsv"
  fi
  touch "$OUT/DONE"
done
n_done=$(ls work/annot/$ID/*/DONE 2>/dev/null | wc -l)
n_exp=$(awk -F'\t' -v i="$ID" 'NR>1 && $1==i && $5=="OK"' work/manifest.tsv | wc -l)
echo -e "STATUS\t$ID\t$n_done/$n_exp assemblies annotated"
