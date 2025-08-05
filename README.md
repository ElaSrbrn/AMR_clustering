# CRE pipeline
# plasmid-outbreak-pipeline

Bioinformatic pipeline for Nanopore plasmid and chromosomal clustering analysis: basecalling → demux → trimming/filtering → assembly → polishing → plasmid & AMR annotation → chromosomal SNP clustering → plasmid clustering.


### `config/config.sh`

```bash
#!/usr/bin/env bash
# -------- USER CONFIG ---------
# basecalling
DORADO_BIN="/path/to/dorado"
CONFIG_FILE="/path/to/dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
INPUT_DIR="/path/to/pod5_dir"   # directory with pod5/fast5
KIT="SQK-RBK114-24"

# project
BASE_DIR="/path/to/project"
OUTDIR="${BASE_DIR}/data"

# barcodes
BCLIST=$(seq -w 01 24)            # or: cat config/barcodes.list

# QC/filtering
QMIN=8
LMIN=200

# assembly/polishing
GENOME_SIZE="5m"
THREADS=8
MODEL="r1041_min_sup_g632"      # medaka model

# reference for chromosomal SNPs
REFERENCE="/path/to/reference.fasta"

# plasmid type for mash matrix (optional example)
PLASMID_TYPE="IncN"
```

### `config/barcodes.list`

```
01
02
...
24
```

### `scripts/utils.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail

log() { echo "[$(date '+%F %T')] $*" >&2; }

check_file() { [[ -s "$1" ]] || { log "Missing file: $1"; return 1; }; }

mkcd() { mkdir -p "$1"; }
```

### `scripts/00_basecall_demux.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
source "$(dirname "$0")/utils.sh"

mkcd "${OUTDIR}"

log "Basecalling with Dorado"
"${DORADO_BIN}" basecaller "${CONFIG_FILE}" "${INPUT_DIR}" \
  --emit-fastq --kit-name "${KIT}" -r --no-trim \
  > "${OUTDIR}/basecalled.fastq"

log "Demultiplexing"
mkcd "${OUTDIR}/demux"
"${DORADO_BIN}" demux \
  --output-dir "${OUTDIR}/demux" \
  --emit-fastq --kit-name "${KIT}" \
  "${OUTDIR}/basecalled.fastq"
```

### `scripts/01_trim_filter.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
source "$(dirname "$0")/utils.sh"

RAW_PREFIX="${OUTDIR}/demux/barcode"   # adjust if dorado naming differs
PC_DIR="${BASE_DIR}/data/porechop"
NF_DIR="${BASE_DIR}/data/filtered_reads"
STATS_DIR="${BASE_DIR}/data/stats"

mkcd "$PC_DIR" "$NF_DIR" "$STATS_DIR"

for bc in ${BCLIST}; do
  in="${RAW_PREFIX}${bc}.fastq"
  [[ -s "$in" ]] || { log "No demux file for $bc"; continue; }

  trim="${PC_DIR}/trimmed_${bc}_passed.fastq"
  out="${NF_DIR}/filtered_barcode${bc}_passed.fastq"

  porechop -i "$in" -o "$trim"
  NanoFilt -q "$QMIN" -l "$LMIN" < "$trim" > "$out"
  rm -f "$trim"
	done

for fq in "${NF_DIR}"/*.fastq; do
  [[ -e "$fq" ]] || continue
  name="${STATS_DIR}/$(basename "${fq%.fastq}").txt"
  seqkit stats -T -a "$fq" > "$name"
done
```

### `scripts/02_assemble_flye.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
source "$(dirname "$0")/utils.sh"

NF_DIR="${BASE_DIR}/data/filtered_reads"
FLY_DIR="${BASE_DIR}/data/unpolished_assemblies"
mkcd "$FLY_DIR"

for bc in ${BCLIST}; do
  reads="${NF_DIR}/filtered_barcode${bc}_passed.fastq"
  outdir="${FLY_DIR}/barcode${bc}"
  [[ -s "$reads" ]] || { log "No reads for $bc"; continue; }

  flye --nano-hq "$reads" \
       --out-dir "$outdir" \
       --genome-size "$GENOME_SIZE" \
       --threads "$THREADS"
done
```

### `scripts/03_polish_medaka.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
source "$(dirname "$0")/utils.sh"

NF_DIR="${BASE_DIR}/data/filtered_reads"
FLY_DIR="${BASE_DIR}/data/unpolished_assemblies"
MIN_DIR="${BASE_DIR}/data/minimap"
MEDAKA_DIR="${BASE_DIR}/data/polished_assemblies"

mkcd "$MIN_DIR" "$MEDAKA_DIR"

for bc in ${BCLIST}; do
  asm="${FLY_DIR}/barcode${bc}/assembly.fasta"
  reads="${NF_DIR}/filtered_barcode${bc}_passed.fastq"
  [[ -s "$asm" ]]   || { log "No assembly for $bc"; continue; }
  [[ -s "$reads" ]] || { log "No reads for $bc"; continue; }

  bam="${MIN_DIR}/aligned_barcode${bc}.bam"
  med_out="${MEDAKA_DIR}/barcode${bc}"
  mkcd "$med_out"

  minimap2 -t "$THREADS" -ax map-ont "$asm" "$reads" \
    | samtools view -b - \
    | samtools sort -@ "$THREADS" -o "$bam" -
  samtools index "$bam"

  medaka_consensus \
    -b "$bam" \
    -d "$asm" \
    -o "$med_out" \
    -t "$THREADS" \
    -m "$MODEL" \
    --bacteria
done
```

### `scripts/04_mob_amr.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
source "$(dirname "$0")/utils.sh"

MEDAKA_DIR="${BASE_DIR}/data/polished_assemblies"
MOB_DIR="${BASE_DIR}/data/plasmid_annotation"
AMR_DIR="${BASE_DIR}/data/amr_annotation"
mkcd "$MOB_DIR" "$AMR_DIR"

for bc in ${BCLIST}; do
  input_fasta="${MEDAKA_DIR}/barcode${bc}/consensus.fasta"
  [[ -s "$input_fasta" ]] || { log "No consensus for $bc"; continue; }

  mob_recon --infile "$input_fasta" --outdir "${MOB_DIR}/barcode${bc}"
  amrfinder -n "$input_fasta" -o "${AMR_DIR}/barcode${bc}.amr.txt" --threads "$THREADS" --plus
done
```

### `scripts/05_snippy_cluster.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
source "$(dirname "$0")/utils.sh"

MEDAKA_DIR="${BASE_DIR}/data/polished_assemblies"
SNIPPY_DIR="${BASE_DIR}/data/snippy"
mkcd "$SNIPPY_DIR"

for SAMPLE in "${MEDAKA_DIR}"/barcode*/consensus.fasta; do
  [[ -s "$SAMPLE" ]] || continue
  [[ "$(basename "$SAMPLE")" == "$(basename "$REFERENCE")" ]] && continue
  BN=$(basename "$(dirname "$SAMPLE")")
  snippy --cpus 4 --ref "$REFERENCE" --ctgs "$SAMPLE" --outdir "${SNIPPY_DIR}/${BN}"
done

snippy-core --ref "$REFERENCE" --prefix snippy_core "${SNIPPY_DIR}"/*
snp-dists snippy_core.full.aln > "${SNIPPY_DIR}/snp_distances.tsv"
```

### `scripts/06_plasmid_cluster.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
source "$(dirname "$0")/utils.sh"

# plasmid.txt: one path per line to plasmid FASTAs from MOB-suite results
mkcd "${BASE_DIR}/data/clustering"
pling align --containment_distance 0.3 --cores 8 --sourmash plasmid.txt pling_out
```

### `scripts/07_mash_distances.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
source "$(dirname "$0")/utils.sh"

output="pairwise_mash_distances_${PLASMID_TYPE}.txt"
> "$output"

for i in {1..24}; do
  for j in {1..24}; do
    if [ "$i" -lt "$j" ]; then
      f1=$(printf "plasmid_%s_barcode%02d.fasta" "$PLASMID_TYPE" "$i")
      f2=$(printf "plasmid_%s_barcode%02d.fasta" "$PLASMID_TYPE" "$j")
      if [[ -f "$f1" && -f "$f2" ]]; then
        mash dist "$f1" "$f2" >> "$output"
      else
        log "Missing file: $f1 or $f2"
      fi
    fi
  done
done
```

