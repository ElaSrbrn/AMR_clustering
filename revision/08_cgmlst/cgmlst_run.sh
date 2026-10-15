#!/usr/bin/env bash
# Submit one SLURM job per scheme for a stage. Run with bash on the login node, not with sbatch.
#   bash cgmlst_run.sh prep [scheme ...]    # default: all schemes in schemes.tsv
#   bash cgmlst_run.sh call [scheme ...]
# Without SLURM: LOCAL=1 bash cgmlst_run.sh prep ecoli   (runs in the foreground)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
stage=${1:?prep|call}; shift
schemes=("$@"); [[ ${#schemes[@]} -gt 0 ]] || mapfile -t schemes < <(awk -F'\t' 'NR>1{print $1}' "$HERE/schemes.tsv")
export HERE
export WORK=${WORK:-$PWD/cgmlst_work}
mkdir -p "$WORK/logs"
for k in "${schemes[@]}"; do
    exp=$(awk -F'\t' -v k="$k" '$1==k{print $3}' "$HERE/schemes.tsv")
    if [[ $exp == CHECK ]]; then echo "skip $k: expected_loci not confirmed in schemes.tsv"; continue; fi
    if [[ $stage == call && ! -d $WORK/$k/inputs ]]; then echo "skip $k: no assemblies in $WORK/$k/inputs"; continue; fi
    if [[ -n ${LOCAL:-} ]]; then bash "$HERE/cgmlst_job.sh" "$stage" "$k"
    else sbatch --parsable --export=ALL -J "cg_${stage}_$k" -o "$WORK/logs/${stage}_$k.%j.out" -e "$WORK/logs/${stage}_$k.%j.err" \
            "$HERE/cgmlst_job.sh" "$stage" "$k" | sed "s/^/submitted $k ($stage): job /"; fi
done
