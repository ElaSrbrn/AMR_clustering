#!/usr/bin/env bash
# run_bart_pipeline.sh  --  read-level MLST + AMR with BART for the nested-depth CRE reads.
# Run with bash on the LOGIN node (setup needs internet; submit calls sbatch itself):
#
#   bash run_bart_pipeline.sh setup     # once: fix env, rebuild MLST schemes, build AMR index, record versions
#   bash run_bart_pipeline.sh submit    # submit one array task per FASTQ (MLST + AMR)
#   bash run_bart_pipeline.sh combine   # after the array has finished: combined + filtered tables
#
# Needs in the same folder: bart_array.sh, bart_job.sh, bart_combine.py, build_bart_amr_index.py

set -euo pipefail

# ============================ CONFIG ==========================================
WORKDIR=/lustre/groups/hpc/urban_lab/projects/plasmid_project/stored_resistant_enterobacterales/new_cre
READ_DIRS=(carblot_pos/10x carblot_pos/20x carblot_neg/10x carblot_neg/20x)   # relative to WORKDIR
BART_ENV=bart
MAX_PARALLEL=20
OUT_PREFIX=bart_results
MIN_REL_DEPTH=0.2             # secondary (truth-independent) AMR filter, see bart_combine.py
# Gold standard for the Cottingham top-n rule: AMRFinder annotation folders of the published
# Flye+Medaka assemblies, highest depth first (first folder containing an isolate is used).
# Set to the real locations of the paper's annotation_<cov> folders.
GOLD_DIRS=(
  carblot_pos/annotation_100x carblot_pos/annotation_80x carblot_pos/annotation_60x
  carblot_neg/annotation_100x carblot_neg/annotation_60x
)
# Same AMRFinderPlus database as the assembly-based annotation
AMRDB=/ictstr01/home/haicu/ela.sauerborn/mamba/envs/amr_env/share/amrfinderplus/data/2026-05-15.1
SCHEME_DIR=$HOME/bart_schemes
# MLST schemes to rebuild: <BART scheme name>|<BIGSdb REST scheme URL>
# (bart-update -p produced empty indexes; names must match 'bart-update -s')
# URLs other than C. freundii are unverified: setup checks every download and stops on failure.
SCHEMES=(
  "Citrobacter_freundii|https://rest.pubmlst.org/db/pubmlst_cfreundii_seqdef/schemes/1"
  "Klebsiella_pneumoniae|https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/schemes/1"
  "Escherichia_coli#1|https://rest.pubmlst.org/db/pubmlst_escherichia_seqdef/schemes/1"
  "Enterobacter_cloacae|https://rest.pubmlst.org/db/pubmlst_ecloacae_seqdef/schemes/1"
  "Klebsiella_oxytoca|https://rest.pubmlst.org/db/pubmlst_koxytoca_seqdef/schemes/1"
)
# ==============================================================================

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$WORKDIR"
bart_db() { micromamba run -n "$BART_ENV" python -c 'import bart,os;print(os.path.dirname(bart.__file__)+"/db")'; }

setup() {
    echo "== 1. Environment fixes (missing dependencies of the Bioconda bart 0.1.2 package)"
    micromamba install -y -n "$BART_ENV" -c conda-forge requests "setuptools<80"
    micromamba run -n "$BART_ENV" python -c "import requests, pkg_resources" \
        || { echo "ERROR: requests/pkg_resources still missing" >&2; exit 1; }
    for t in kma refseq_masher mash; do
        micromamba run -n "$BART_ENV" which "$t" >/dev/null || { echo "ERROR: $t missing in env $BART_ENV" >&2; exit 1; }
    done
    DB=$(bart_db); mkdir -p "$SCHEME_DIR"

    echo "== 2. Rebuild MLST schemes from BIGSdb REST"
    for entry in "${SCHEMES[@]}"; do
        name="${entry%%|*}"; url="${entry#*|}"
        echo "-- $name"
        curl -sf "$url/profiles_csv" > "$SCHEME_DIR/$name.tab" || { echo "ERROR: download failed: $url/profiles_csv" >&2; exit 1; }
        head -1 "$SCHEME_DIR/$name.tab" | grep -q "^ST"$'\t' || { echo "ERROR: $name profile table does not start with ST: $(head -c 120 "$SCHEME_DIR/$name.tab")" >&2; exit 1; }
        curl -sf "$url" | micromamba run -n "$BART_ENV" python -c 'import sys,json;[print(l) for l in json.load(sys.stdin)["loci"]]' \
            | while read -r l; do curl -sf "$l/alleles_fasta"; done > "$SCHEME_DIR/$name.fna"
        n_loci=$(head -1 "$SCHEME_DIR/$name.tab" | awk -F'\t' '{n=0; for(i=2;i<=NF;i++) if($i!="clonal_complex" && $i!="species" && $i!="CC" && $i!="Lineage") n++; print n}')
        n_al=$(grep -c ">" "$SCHEME_DIR/$name.fna" || true)
        (( n_al > n_loci )) || { echo "ERROR: $name has only $n_al alleles" >&2; exit 1; }
        echo "   $n_loci loci, $n_al alleles, $(($(wc -l < "$SCHEME_DIR/$name.tab")-1)) STs"
        ( cd "$SCHEME_DIR" && micromamba run -n "$BART_ENV" bart-update -a "$name.fna" "$name.tab" >/dev/null 2>&1 )
        [[ -s "$DB/indexes/$name.name" ]] || { echo "ERROR: KMA index for $name not built" >&2; exit 1; }
    done

    echo "== 3. AMR index from AMRFinderPlus database $(basename "$AMRDB")"
    amr_fa="$SCHEME_DIR/amr_$(basename "$AMRDB").fa"
    micromamba run -n "$BART_ENV" python3 "$HERE/build_bart_amr_index.py" "$AMRDB" "$amr_fa" \
        --symbol-field "${AMR_SYMBOL_FIELD:-6}" --desc-field "${AMR_DESC_FIELD:-8}"
    micromamba run -n "$BART_ENV" kma index -i "$amr_fa" -o "$DB/indexes/amr"
    [[ -s "$DB/indexes/amr.name" ]] || { echo "ERROR: AMR index not built" >&2; exit 1; }

    echo "== 4. Record versions"
    {
        echo "date: $(date -Iseconds)"
        echo "bart: $(micromamba run -n "$BART_ENV" bart -h 2>&1 | grep -m1 -oi 'bart [0-9.]*' || echo 0.1.2)"
        echo "kma: $(micromamba run -n "$BART_ENV" kma -v 2>&1 | head -1)"
        echo "amrfinder_db: $AMRDB"
        for entry in "${SCHEMES[@]}"; do echo "mlst_scheme: ${entry%%|*} <- ${entry#*|} ($(date +%F))"; done
        echo "amr_index_templates: $(wc -l < "$DB/indexes/amr.name")"
    } > "$WORKDIR/bart_setup_info.txt"
    micromamba env export -n "$BART_ENV" > "$WORKDIR/bart_env_export.yml"
    echo "Setup done: $WORKDIR/bart_setup_info.txt, $WORKDIR/bart_env_export.yml"
}

submit() {
    [[ -f "$WORKDIR/bart_setup_info.txt" ]] || { echo "ERROR: run 'setup' first" >&2; exit 1; }
    for f in bart_array.sh bart_job.sh; do [[ -f "$HERE/$f" ]] || { echo "ERROR: missing $HERE/$f" >&2; exit 1; }; done
    BART_ENV="$BART_ENV" bash "$HERE/bart_array.sh" -p "$MAX_PARALLEL" "${READ_DIRS[@]}"
}

combine() {
    gold=()
    for g in "${GOLD_DIRS[@]}"; do
        if [[ -d "$g" ]]; then gold+=("$g"); else echo "NOTE: gold dir not found, skipped: $g" >&2; fi
    done
    (( ${#gold[@]} > 0 )) || { echo "ERROR: none of GOLD_DIRS exists; set them in the CONFIG block" >&2; exit 1; }
    python3 "$HERE/bart_combine.py" "$OUT_PREFIX" "${READ_DIRS[@]}" \
        --catalog "$AMRDB/ReferenceGeneCatalog.txt" --min-rel-depth "$MIN_REL_DEPTH" \
        --gold-dirs "${gold[@]}"
}

case "${1:-}" in
    setup) setup ;;
    submit) submit ;;
    combine) combine ;;
    *) echo "Usage: bash $0 setup|submit|combine" >&2; exit 1 ;;
esac
