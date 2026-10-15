#!/bin/bash
#SBATCH --job-name=pcs_gate
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:15:00
#SBATCH --exclude=cpusrv08
#SBATCH --output=logs/01b_gate_%j.out
# Finish step 01 unattended, then submit steps 02-04. Usage: bash scripts/01b_finish_annotation.sh
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
ROUND=${ROUND:-1}
MAX_ROUNDS=${MAX_ROUNDS:-3}
source config.sh

running=$(squeue --me -h -n pcs_annotate -o %A | sort -u | paste -sd:)
if [ -n "$running" ]; then
  gid=$(sbatch --parsable --dependency=afterany:"$running" --export=ALL,ROUND="$ROUND",MAX_ROUNDS="$MAX_ROUNDS" \
        scripts/01b_finish_annotation.sh)
  echo "annotation jobs $running still active; check queued as job $gid"
  exit 0
fi

missing=()
while read -r ID; do
  [ -z "$ID" ] && continue
  exp=$(awk -F'\t' -v i="$ID" 'NR>1 && $1==i && $5=="OK"' work/manifest.tsv | wc -l)
  got=$(ls work/annot/$ID/*/DONE 2>/dev/null | wc -l)
  if [ "$exp" -ne "$got" ]; then
    idx=$(grep -nx "$ID" "$IDS_FILE" | cut -d: -f1)
    missing+=("$idx")
    echo "round $ROUND: $ID (index $idx) $got/$exp annotated"
  fi
done < "$IDS_FILE"

if [ ${#missing[@]} -eq 0 ]; then
  echo "round $ROUND: all assemblies annotated"
  jid=$(sbatch --parsable scripts/02_03_collect_compare.sbatch)
  echo "submitted steps 02-04 as job $jid  (log: logs/02_03_compare_${jid}.out)"
  exit 0
fi

if [ "$ROUND" -gt "$MAX_ROUNDS" ]; then
  echo "still incomplete after $MAX_ROUNDS rounds; check these logs:"
  for i in "${missing[@]}"; do ls -t logs/01_annotate_*_"$i".out 2>/dev/null | head -1; done
  exit 1
fi

list=$(IFS=,; echo "${missing[*]}")
aid=$(sbatch --parsable --array="${list}%3" scripts/01_annotate.sbatch)
gid=$(sbatch --parsable --dependency=afterany:"$aid" --export=ALL,ROUND=$((ROUND + 1)),MAX_ROUNDS="$MAX_ROUNDS" \
      scripts/01b_finish_annotation.sh)
echo "round $ROUND: resubmitted indices $list as array $aid; next check is job $gid"
