
#!/bin/bash
#SBATCH --job-name=downsample_annot_cpu
#SBATCH --partition=cpu_p
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --mail-user=ela.sauerborn@helmholtz-munich.de
#SBATCH --mail-type=ALL
#SBATCH --nice=10000


# -------- settings --------
MERGED_FASTQ="merged_fastq"         # dir with barcode0102_L1000_Q9.fastq.gz
OUT="assemblies/by_coverage"        # output root
GENOME_BP=5000000                   # ~5 Mbp; adjust per organism
DEPTHS="$(seq 10 10 100)"           # 10x,20x,...,100x
THREADS=16
SEED=42
POLISH=1                            # 0 to skip Medaka polishing

# envs
ASM_ENV="assembly"      # flye minimap2 medaka samtools seqtk
MOB_ENV="mobsuite_env"  # mob_suite=3.1.9
AMR_ENV="amr_env"       # ncbi-amrfinderplus=4.0.23
MLST_ENV="mlst_env"     # mlst

mkdir -p "$OUT"
# Optional, safe to repeat (skip if offline)
# micromamba run -n "$AMR_ENV" amrfinder --update || true

process_one() {
  local INPUT_FASTQ="$1"
  local tag base cov covdir subset draft polished asm_for_anno
  tag="$(basename "$INPUT_FASTQ")"
  tag="${tag%_L1000_Q9.fastq.gz}"  # e.g., barcode0102

  local SUBDIR="$OUT/$tag"
  local SHUF_DIR="$SUBDIR/shuffle"
  local SUBSETS_DIR="$SUBDIR/subsets"
  mkdir -p "$SHUF_DIR" "$SUBSETS_DIR"

  # Shuffle once (seqtk reads gz)
  local SHUF_FQ="$SHUF_DIR/${tag}.shuffled.fastq"
  if [[ ! -s "$SHUF_FQ" ]]; then
    echo "[$tag] shuffling reads…"
    micromamba run -n "$ASM_ENV" seqtk sample -s "$SEED" "$INPUT_FASTQ" 1 > "$SHUF_FQ"
  fi

  # Target bases per depth and read-counts needed
  local TARGETS_TSV="$SHUF_DIR/targets.tsv"; : > "$TARGETS_TSV"
  for d in $DEPTHS; do echo -e "${d}\t$((GENOME_BP * d))" >> "$TARGETS_TSV"; done

  local COUNTS_TSV="$SHUF_DIR/depth_read_counts.tsv"
  awk -v OFS='\t' '
    NR==FNR { tgt[$1]=$2; depths[n++]=$1; next }
    (NR%4)==2 { cum+=length($0); reads+=1;
      while (i<n && cum>=tgt[depths[i]]) { print depths[i], reads; i++ } }
    END { if (i<n) { printf("[WARN] %s insufficient for:", ARGV[2]) > "/dev/stderr";
      for (; i<n; i++) printf(" %sx", depths[i]) > "/dev/stderr"; printf("\n") > "/dev/stderr" } }
  ' "$TARGETS_TSV" "$SHUF_FQ" > "$COUNTS_TSV"

  # Materialize subsets
  while IFS=$'\t' read -r depth reads_needed; do
    subset="$SUBSETS_DIR/${tag}.depth_${depth}x.fastq"
    [[ -s "$subset" ]] && continue
    echo "[$tag] depth ${depth}× → ${reads_needed} reads"
    head -n $((reads_needed * 4)) "$SHUF_FQ" > "$subset"
  done < "$COUNTS_TSV"

  # Assemble + (optional) polish + annotate per subset
  while read -r subset; do
    [[ -s "$subset" ]] || continue
    base="$(basename "$subset" .fastq)"      # barcode0102.depth_40x
    cov="${base##*.depth_}"; cov="${cov%x}" # "40"
    covdir="$SUBDIR/cov${cov}"
    mkdir -p "$covdir"

    draft="$covdir/assembly_cov${cov}.fasta"
    if [[ ! -s "$draft" ]]; then
      echo "[$tag] cov ${cov}: Flye"
      micromamba run -n "$ASM_ENV" flye \
        --nano-hq "$subset" \
        --genome-size "${GENOME_BP}" \
        --threads "$THREADS" \
        --plasmids \
        --out-dir "$covdir"
      cp "$covdir/assembly.fasta" "$draft"
    fi

    asm_for_anno="$draft"
    if [[ "$POLISH" -eq 1 ]]; then
      polished="$covdir/polished_cov${cov}.fasta"
      if [[ ! -s "$polished" ]]; then
        echo "[$tag] cov ${cov}: Medaka --bacteria (CPU)"
        micromamba run -n "$ASM_ENV" minimap2 -x map-ont -t "$THREADS" "$draft" "$subset" \
          | micromamba run -n "$ASM_ENV" samtools sort -@ "$THREADS" -o "$covdir/reads.sorted.bam"
        micromamba run -n "$ASM_ENV" samtools index "$covdir/reads.sorted.bam"

        micromamba run -n "$ASM_ENV" medaka consensus \
          -i "$subset" -d "$draft" -o "$covdir/medaka" -t "$THREADS" --bacteria
        cp "$covdir/medaka/consensus.fasta" "$polished"
        rm -rf "$covdir/medaka" "$covdir/reads.sorted.bam" "$covdir/reads.sorted.bam.bai"
      fi
      asm_for_anno="$polished"
    fi

    # MOB-suite
    if [[ ! -d "$covdir/mobsuite" ]]; then
      echo "[$tag] cov ${cov}: MOB-suite"
      micromamba run -n "$MOB_ENV" mob_recon -i "$asm_for_anno" -o "$covdir/mobsuite" \
        --run_typer --run_mobtyper
    fi

    # AMRFinderPlus
    if [[ ! -s "$covdir/amrfinder_cov${cov}.tsv" ]]; then
      echo "[$tag] cov ${cov}: AMRFinderPlus"
      micromamba run -n "$AMR_ENV" amrfinder \
        -n "$asm_for_anno" -O bacteria --threads "$THREADS" \
        -o "$covdir/amrfinder_cov${cov}.tsv"
    fi

    # MLST (polished if present, else draft)
    if [[ ! -s "$covdir/mlst_cov${cov}.tsv" ]]; then
      micromamba run -n "$MLST_ENV" mlst "$asm_for_anno" > "$covdir/mlst_cov${cov}.tsv" || true
    fi
  done < <(find "$SUBSETS_DIR" -name "${tag}.depth_*x.fastq" -type f | sort -V)

  echo "[$tag] done → $SUBDIR"
}

shopt -s nullglob
for fq in "$MERGED_FASTQ"/barcode*_L1000_Q9.fastq.gz; do
  process_one "$fq"
done

