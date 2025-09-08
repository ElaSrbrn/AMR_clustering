READS_DIR="/path/to/raw_pod5_or_fast5"   # directory with raw ONT reads
KIT_NAME="SQK-RBK114-96"                 # your kit
MODELS_DIR="/lustre/groups/hpc/urban_lab/tools/dorado-0.9.1-linux-x64/bin"
OUTDIR="basecalled_v5_mod"

mkdir -p "$OUTDIR"
./dorado basecaller \
  --kit-name "$KIT_NAME" \
  --models-directory "$MODELS_DIR" \
  --trim \
  --emit-bam \
  sup 6mA 4mC_5mC \
  -r "$READS_DIR" > "$OUTDIR/all.bam"

mkdir -p "$OUTDIR/demux_bam"
./dorado demux \
  --kit-name "$KIT_NAME" \
  --output-dir "$OUTDIR/demux_bam" \
  "$OUTDIR/all.bam"

mkdir -p "$OUTDIR/demux_fastq"
for b in "$OUTDIR"/demux_bam/*.bam; do
  base=$(basename "${b%.bam}")
  samtools fastq "$b" | gzip > "$OUTDIR/demux_fastq/${base}.fastq.gz"
done
