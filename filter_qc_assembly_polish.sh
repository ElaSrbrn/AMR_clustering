#!/bin/bash
#SBATCH -p gpu_p
#SBATCH -q gpu
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH --qos gpu_normal
#SBATCH --partition=gpu_p
#SBATCH --nice=10000
#SBATCH --mail-user=ela.sauerborn@helmholtz-munich.de
#SBATCH --mail-type=ALL
#SBATCH --gres=gpu:2
#SBATCH --job-name=assembly
#SBATCH -c 2

# ==== INPUTS (from your basecalling script) ====
OUTDIR="basecalled_v5_mod"          # where demux_fastq/ lives
DEMUX_FASTQ="$OUTDIR/demux_fastq"

# ==== ENV / TOOLS ====
ASSEMBLY_ENV="assembly_env"         # chopper seqkit flye minimap2 medaka samtools
MOB_ENV="mobsuite_env"              # mob_suite 3.1.9
AMR_ENV="amr_env"                   # ncbi-amrfinderplus 4.0.23
MLST_ENV="mlst_env"                 # mlst
THREADS=16
GENOME_SIZE="5m"                    # ~ genome size (e.g., 5m for bacteria)

# (optional) MLST database path; leave empty to use default bundled DB
MLST_DB=""

# ==== OUTPUTS ====
FILTERED_FASTQ="$OUTDIR/filtered_fastq"
MERGED_FASTQ="$OUTDIR/merged_fastq"
ASM_UNPOLISHED="assemblies/unpolished"
ASM_POLISHED="assemblies/polished"

mkdir -p "$FILTERED_FASTQ" "$MERGED_FASTQ" "$ASM_UNPOLISHED" "$ASM_POLISHED"

echo "==> 1) Filter each barcode (len>=1000, Q>=9)"
shopt -s nullglob
for fq in "$DEMUX_FASTQ"/*.fastq.gz; do
  base=$(basename "${fq%.fastq.gz}")                      # barcodeXX
  out="$FILTERED_FASTQ/${base}_L1000_Q9.fastq.gz"
  [[ -s "$out" ]] && continue
  micromamba run -n "$ASSEMBLY_ENV" chopper -t "$THREADS" -l 1000 -q 9 -i "$fq" \
    | gzip > "$out"
done

echo "==> 2) Concatenate pairs (01+02→0102 … 95+96→9596) and deduplicate read IDs"
for ((i=1; i<=96; i+=2)); do
  b1=$(printf "barcode%02d_L1000_Q9.fastq.gz" "$i")
  b2=$(printf "barcode%02d_L1000_Q9.fastq.gz" "$((i+1))")
  out=$(printf "barcode%02d%02d_L1000_Q9.fastq.gz" "$i" "$((i+1))")
  f1="$FILTERED_FASTQ/$b1"
  f2="$FILTERED_FASTQ/$b2"
  fout="$MERGED_FASTQ/$out"

  [[ -f "$f1" && -f "$f2" ]] || { echo "skip pair $b1 + $b2"; continue; }
  [[ -s "$fout" ]] && continue

  if micromamba run -n "$ASSEMBLY_ENV" which seqkit >/dev/null 2>&1; then
    micromamba run -n "$ASSEMBLY_ENV" bash -lc \
      "zcat \"$f1\" \"$f2\" | seqkit rmdup -n | gzip > \"$fout\""
  else
    zcat "$f1" "$f2" | awk '
      NR%4==1 { id=$1; sub(/^@/,"",id); skip=seen[id]++ }
      { buf[NR%4]=$0; if (NR%4==0 && !skip) print buf[1] ORS buf[2] ORS buf[3] ORS buf[0] }
    ' | gzip > "$fout"
  fi
done

echo "==> 3) Assemble (Flye --nano-hq) and 4) Polish (Medaka --bacteria)"
for fq in "$MERGED_FASTQ"/barcode*_L1000_Q9.fastq.gz; do
  [[ -e "$fq" ]] || continue
  sample=$(basename "${fq%_L1000_Q9.fastq.gz}")          # e.g., barcode0102
  draft="$ASM_UNPOLISHED/${sample}.fasta"
  polished="$ASM_POLISHED/${sample}.fasta"

  if [[ ! -s "$draft" ]]; then
    tmp_flye=$(mktemp -d)
    micromamba run -n "$ASSEMBLY_ENV" flye \
      --nano-hq "$fq" \
      --out-dir "$tmp_flye" \
      --threads "$THREADS" \
      --genome-size "$GENOME_SIZE" \
      --plasmids \
      --asm-coverage 60
    cp "$tmp_flye/assembly.fasta" "$draft"
    rm -rf "$tmp_flye"
  fi

  if [[ ! -s "$polished" ]]; then
    micromamba run -n "$ASSEMBLY_ENV" minimap2 -x map-ont -t "$THREADS" "$draft" "$fq" \
      | micromamba run -n "$ASSEMBLY_ENV" samtools sort -@ "$THREADS" -o "${draft%.fasta}.reads.bam"
    micromamba run -n "$ASSEMBLY_ENV" samtools index "${draft%.fasta}.reads.bam"

    tmp_medaka=$(mktemp -d)
    micromamba run -n "$ASSEMBLY_ENV" medaka consensus \
      -i "$fq" \
      -d "$draft" \
      -o "$tmp_medaka" \
      -t "$THREADS" \
      --bacteria
    cp "$tmp_medaka/consensus.fasta" "$polished"
    rm -rf "$tmp_medaka" "${draft%.fasta}.reads.bam" "${draft%.fasta}.reads.bam.bai"
  fi

  echo "→ ${sample}: unpolished=$draft  polished=$polished"

  # ==== 5) Annotations on POLISHED assemblies ====

  # MOB-suite (needs a directory)
  mob_dir="$ASM_POLISHED/${sample}_mobsuite"
  if [[ ! -d "$mob_dir" ]]; then
    micromamba run -n "$MOB_ENV" mob_recon \
      -i "$polished" \
      -o "$mob_dir" \
      --run_typer --run_mobtyper
  fi

  # AMRFinderPlus (single TSV)
  amr_tsv="$ASM_POLISHED/${sample}.amrfinder.tsv"
  if [[ ! -s "$amr_tsv" ]]; then
    micromamba run -n "$AMR_ENV" amrfinder \
      -n "$polished" \
      -O bacteria \
      --threads "$THREADS" \
      -o "$amr_tsv"
  fi

  # MLST (single TSV) – use MLST_DB if set
  mlst_tsv="$ASM_POLISHED/${sample}.mlst.tsv"
  if [[ ! -s "$mlst_tsv" ]]; then
    if [[ -n "${MLST_DB}" ]]; then
      micromamba run -n "$MLST_ENV" mlst --datadir "$MLST_DB" "$polished" > "$mlst_tsv"
    else
      micromamba run -n "$MLST_ENV" mlst "$polished" > "$mlst_tsv"
    fi
  fi
done

echo "All done."
