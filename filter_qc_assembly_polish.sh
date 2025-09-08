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
ASSEMBLY_ENV="assembly_env"         # has: chopper seqkit flye minimap2 medaka samtools
THREADS=16
GENOME_SIZE="5m"                    # ~ genome size (e.g., 5m for bacteria)

# ==== OUTPUTS ====
FILTERED_FASTQ="$OUTDIR/filtered_fastq"
MERGED_FASTQ="$OUTDIR/merged_fastq"
ASM_UNPOLISHED="assemblies/unpolished"
ASM_POLISHED="assemblies/polished"

mkdir -p "$FILTERED_FASTQ" "$MERGED_FASTQ" "$ASM_UNPOLISHED" "$ASM_POLISHED"

echo "==> 1) Filter each barcode (len>=1000, Q>=9)"
for fq in "$DEMUX_FASTQ"/*.fastq.gz; do
  base=$(basename "${fq%.fastq.gz}")                      # barcodeXX
  micromamba run -n "$ASSEMBLY_ENV" chopper -t "$THREADS" -l 1000 -q 9 -i "$fq" \
    | gzip > "$FILTERED_FASTQ/${base}_L1000_Q9.fastq.gz"
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

  if micromamba run -n "$ASSEMBLY_ENV" which seqkit >/dev/null 2>&1; then
    micromamba run -n "$ASSEMBLY_ENV" bash -lc \
      "zcat \"$f1\" \"$f2\" | seqkit rmdup -n | gzip > \"$fout\""
  else
    # fallback dedup by read ID
    zcat "$f1" "$f2" | awk '
      NR%4==1 { id=$1; sub(/^@/,"",id); skip=seen[id]++ }
      { buf[NR%4]=$0; if (NR%4==0 && !skip) print buf[1] ORS buf[2] ORS buf[3] ORS buf[0] }
    ' | gzip > "$fout"
  fi
done

echo "==> 3) Assemble (Flye --nano-hq) and 4) Polish (Medaka --bacteria) — flat outputs"
for fq in "$MERGED_FASTQ"/barcode*_L1000_Q9.fastq.gz; do
  [[ -e "$fq" ]] || continue
  sample=$(basename "${fq%_L1000_Q9.fastq.gz}")          # e.g., barcode0102
  draft="$ASM_UNPOLISHED/${sample}.fasta"
  polished="$ASM_POLISHED/${sample}.fasta"

  # Flye requires an output dir; use a temp and collapse to a single FASTA
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

  # Map reads to draft and polish with Medaka (also temp dir → single FASTA)
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

  echo "→ ${sample}: unpolished=$draft  polished=$polished"
done

echo "Done."
