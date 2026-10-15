# ---------------------------------------------------------------------------
# Plasmid contig stability pipeline: configuration
# Source this file in every step:  source config.sh
# Lines marked PLACEHOLDER must be filled in before running.
# ---------------------------------------------------------------------------

# Locations
BASE=/ictstr01/groups/hpc/urban_lab/projects/plasmid_project
PCS=$BASE/plasmid_contig_stability                  # this folder
NEWCRE=$BASE/stored_resistant_enterobacterales/new_cre
ISOLATE_DIRS=$BASE/revision_scripts/isolate_dirs_S1.txt   # Autocycler isolate folders (60x)
IDS_FILE=$PCS/ids_s1a_60x.txt                       # carbapenemase-positive isolates with a 60x reference (98)
HYB_OUT=$BASE/hybracter_depth_out_medaka            # Hybracter (Medaka-polished) output
DEPTHS="10x 20x 40x 60x"

# Flye assemblies (Medaka-polished), one FASTA per isolate and depth.
# Shell glob with {ID} and {DEPTH}; must match exactly one file per isolate/depth (6901: see overrides.tsv).
# Circularity is read from the matching flye_<depth>/<name>_flye/assembly_info.txt.
FLYE_GLOB='/ictstr01/groups/hpc/urban_lab/projects/plasmid_project/stored_resistant_enterobacterales/new_cre/carblot_*/polished_flye_{DEPTH}/isolate_{ID}_*_medaka/consensus.fasta'

# Environments
ENV_MAP=autocycler_0.8.0      # minimap2, samtools, python3
ENV_AMR=amr_env               # AMRFinderPlus (also provides BLAST+)
ENV_MOB=mobsuite_env          # MOB-typer
ENV_BAKTA=bakta_clean         # Bakta v1.11.4
ENV_BLAST=amr_env             # blastn for ISfinder and self-repeats

# Databases
BAKTA_DB=/ictstr01/groups/hpc/urban_lab/backup/bakta_db_v6/db   # as in the EMBL project
ISFINDER_DB=/ictstr01/home/haicu/ela.sauerborn/EMBL/ggCallaroo/ISfinder_db   # as in the EMBL project

# ISfinder hit filters. PLACEHOLDER values: replace with the thresholds used in the EMBL project.
IS_MIN_ID=95                  # % identity
IS_MIN_COV=80                 # % of the IS element covered
REPEAT_MIN_LEN=500            # self-repeats: minimum length (bp)
REPEAT_MIN_ID=90              # self-repeats: minimum identity (%)

# Stability criteria (from the Flye vs Autocycler comparison document)
LEN_TOL=0.001                 # |length difference| / reference length
MIN_AF=0.99                   # fraction of the reference aligned
MIN_ID=0.999                  # identity across aligned blocks
MIN_BLOCK=1000                # ignore alignment blocks shorter than this (bp)
BREAK_GAP=1000                # discontinuity (bp) that counts as a breakpoint
CHROM_MIN=1000000             # contigs at least this long are treated as chromosomal
CHIMERA_UNALIGNED=0.05        # primary contig with more than this fraction unaligned = chimeric
FRAGMENT_MIN_BP=2000          # other contigs with at least this much aligned sequence count as fragments
SMALL_PLASMID=20000           # plasmids below this length are proposed for synteny by default
IS_WINDOW=1000                # breakpoint within this distance of an IS element or repeat = IS/repeat-associated

# Python from the mapping environment (stdlib only)
PYBIN="${MAMBA_ROOT_PREFIX:-$HOME/mamba}/envs/$ENV_MAP/bin/python3"; if [ -x "$PYBIN" ]; then PY="$PYBIN"; else PY="micromamba run -n $ENV_MAP python3"; fi

# Node-local temporary directory (the default /tmp filled up on some nodes)
make_jobtmp () {
  JOBTMP=/localscratch/$USER/pcs_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}
  if ! { mkdir -p "$JOBTMP" 2>/dev/null && [ -d "$JOBTMP" ] && [ -w "$JOBTMP" ]; }; then
    echo "WARN: $JOBTMP not usable on $(hostname); falling back to Lustre"
    JOBTMP=$PCS/tmp_jobs/pcs_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}
    mkdir -p "$JOBTMP" || { echo "FATAL: no usable temporary directory"; exit 1; }
  fi
  export TMPDIR=$JOBTMP TMP=$JOBTMP TEMP=$JOBTMP
  trap 'rm -rf "$JOBTMP"' EXIT
}
