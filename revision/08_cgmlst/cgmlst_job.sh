#!/usr/bin/env bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH -t 24:00:00
#SBATCH --nice=10000
#SBATCH --job-name=cgmlst
# One stage for one scheme. Submitted by cgmlst_run.sh; can also be run directly with bash.
#   bash cgmlst_job.sh prep <scheme>   unzip scheme, Prodigal training file, PrepExternalSchema
#   bash cgmlst_job.sh call <scheme>   AlleleCall on all assemblies (all depths) of that scheme
set -eo pipefail
stage=${1:?prep|call}; k=${2:?scheme key}
WORK=${WORK:-cgmlst_work}; RAW=${RAW:-$HOME/cgmlst/raw}; ENV=${CHEWIE_ENV:-chewie}
CPU=${SLURM_CPUS_PER_TASK:-8}
eval "$(micromamba shell hook --shell=bash)"; micromamba activate "$ENV"; set -u
SCRIPTS=${HERE:-$(cd "$(dirname "$0")" && pwd)}
row=$(awk -F'\t' -v k="$k" '$1==k' "$SCRIPTS/schemes.tsv")
[[ -n "$row" ]] || { echo "ERROR: $k not in schemes.tsv" >&2; exit 1; }
IFS=$'\t' read -r _ sid expected _ <<< "$row"
d="$WORK/$k"; mkdir -p "$d"

if [[ $stage == prep ]]; then
    zip=""; for z in "$RAW/$sid.zip" "$RAW/$k.zip"; do [[ -f $z ]] && { zip=$z; break; }; done
    [[ -n $zip ]] || { echo "ERROR: no zip for $k ($RAW/$sid.zip)" >&2; exit 1; }
    rm -rf "$d/raw"; mkdir -p "$d/raw"; unzip -q -o "$zip" -d "$d/raw"
    n=$(find "$d/raw" -type f \( -name '*.fasta' -o -name '*.fa' -o -name '*.fna' \) | wc -l)
    echo "$k: scheme $sid, $n locus files (expected $expected)"
    if [[ $expected != "$n" ]]; then
        echo "ERROR: locus count $n does not match expected_loci ($expected) in schemes.tsv; check the scheme page" >&2; exit 1
    fi
    # chewBBACA expects one FASTA per locus in one folder
    if [[ -z $(find "$d/raw" -maxdepth 1 -type f -name '*.fasta' | head -1) ]]; then
        mkdir -p "$d/raw_flat"; find "$d/raw" -type f \( -name '*.fasta' -o -name '*.fa' -o -name '*.fna' \) -exec ln -sf {} "$d/raw_flat/" \;
        src="$d/raw_flat"
    else src="$d/raw"; fi
    tg=$(awk -F'\t' -v k="$k" '$1==k{print $3}' "$WORK/training.tsv")
    [[ -f $tg ]] || { echo "ERROR: no training genome for $k in $WORK/training.tsv (run cgmlst_manifest.py)" >&2; exit 1; }
    prodigal -i "$tg" -t "$d/$k.trn" -p single > /dev/null
    rm -rf "$d/schema"
    chewBBACA.py PrepExternalSchema -g "$src" -o "$d/schema" --ptf "$d/$k.trn" --cpu "$CPU"
    { echo "scheme_key: $k"; echo "cgmlst_org_id: $sid"; echo "zip: $zip ($(date -r "$zip" +%F))"
      echo "loci_in_zip: $n"; echo "training_genome: $tg"; echo "prepared: $(date -Iseconds)"
      echo "chewBBACA: $(chewBBACA.py --version 2>&1 | tail -1)"; echo "prodigal: $(prodigal -v 2>&1 | grep -i prodigal | head -1)"
    } > "$d/prep_info.txt"
    echo "done: $d/schema"

elif [[ $stage == call ]]; then
    [[ -d $d/schema ]] || { echo "ERROR: $d/schema missing; run prep first" >&2; exit 1; }
    n_in=$(ls "$d/inputs" | wc -l); echo "$k: $n_in assemblies"
    # work on a copy: AlleleCall adds inferred alleles to the schema it is given
    rm -rf "$d/schema_run" "$d/allelecall"; cp -r "$d/schema" "$d/schema_run"
    chewBBACA.py AlleleCall -i "$d/inputs" -g "$d/schema_run" -o "$d/allelecall" --cpu "$CPU" --output-novel --force-continue
    echo "done: $d/allelecall/results_alleles.tsv"
else
    echo "ERROR: stage must be prep or call" >&2; exit 1
fi
