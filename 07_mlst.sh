#!/usr/bin/env bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=32G
#SBATCH -t 08:00:00
#SBATCH --nice=10000

set -euo pipefail
shopt -s nullglob

INPUT_GLOB="${INPUT_GLOB:-polished_fasta_*}"
OUTROOT="${OUTROOT:-taxonomy_mlst_results}"
MLST_ENV="${MLST_ENV:-mlst_env}"

mkdir -p "${OUTROOT}"

OVERALL_TSV="${OUTROOT}/all_coverages_mlst.tsv"
OVERALL_CSV="${OUTROOT}/all_coverages_mlst.csv"

printf "coverage\tisolate\tassembly_path\tmlst_scheme\tmlst_st\tmlst_alleles\n" > "${OVERALL_TSV}"

have_any_dir=0

for INDIR in ${INPUT_GLOB}; do
  [[ -d "${INDIR}" ]] || continue
  have_any_dir=1

  base="$(basename "${INDIR}")"
  cov="$(printf '%s\n' "${base}" | grep -Eo '[0-9]+x' | head -n1 || true)"
  if [[ -z "${cov}" ]]; then
    cov="${base}"
  fi

  OUTDIR="${OUTROOT}/${cov}"
  MLST_DIR="${OUTDIR}/mlst"
  mkdir -p "${MLST_DIR}"

  MANIFEST="${OUTDIR}/assemblies_manifest.tsv"
  FOFN="${OUTDIR}/assemblies.fofn"
  MLST_TSV="${MLST_DIR}/mlst.tsv"
  MLST_NORM="${MLST_DIR}/mlst_normalized.tsv"
  COMBINED_TSV="${OUTDIR}/${cov}_mlst.tsv"
  COMBINED_CSV="${OUTDIR}/${cov}_mlst.csv"

  : > "${MANIFEST}"
  : > "${FOFN}"

  while IFS= read -r -d '' f; do
    bn="$(basename "${f}")"
    isolate="${bn}"
    isolate="${isolate%.fasta.gz}"
    isolate="${isolate%.fa.gz}"
    isolate="${isolate%.fna.gz}"
    isolate="${isolate%.fasta}"
    isolate="${isolate%.fa}"
    isolate="${isolate%.fna}"

    printf "%s\t%s\t%s\n" "${isolate}" "${f}" "${cov}" >> "${MANIFEST}"
    printf "%s\n" "${f}" >> "${FOFN}"
  done < <(
    find -L "${INDIR}" -maxdepth 1 \( -type f -o -type l \) \
      \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \
         -o -name "*.fasta.gz" -o -name "*.fa.gz" -o -name "*.fna.gz" \) \
      -print0 | sort -z
  )

  if [[ ! -s "${MANIFEST}" ]]; then
    echo "[WARN] No FASTA assemblies found in ${INDIR}; skipping"
    continue
  fi

  echo "[INFO] Processing ${INDIR} -> coverage ${cov}"

  mapfile -t assemblies < "${FOFN}"
  micromamba run -n "${MLST_ENV}" mlst "${assemblies[@]}" > "${MLST_TSV}"

  awk -F'\t' '
    BEGIN {
      OFS="\t"
      print "isolate","mlst_scheme","mlst_st","mlst_alleles"
    }
    {
      file=$1
      scheme=(NF>=2 ? $2 : "")
      st=(NF>=3 ? $3 : "")
      alleles=""
      if (NF>=4) {
        alleles=$4
        for (i=5; i<=NF; i++) alleles=alleles ";" $i
      }

      n=file
      sub(/^.*\//, "", n)
      sub(/\.fasta\.gz$/, "", n)
      sub(/\.fa\.gz$/, "", n)
      sub(/\.fna\.gz$/, "", n)
      sub(/\.fasta$/, "", n)
      sub(/\.fa$/, "", n)
      sub(/\.fna$/, "", n)

      print n, scheme, st, alleles
    }
  ' "${MLST_TSV}" > "${MLST_NORM}"

  awk -F'\t' '
    BEGIN {
      OFS="\t"
      print "coverage","isolate","assembly_path","mlst_scheme","mlst_st","mlst_alleles"
    }
    FNR==NR {
      if (FNR>1) {
        mlst_scheme[$1]=$2
        mlst_st[$1]=$3
        mlst_alleles[$1]=$4
      }
      next
    }
    FNR>1 {
      isolate=$1
      path=$2
      cov=$3
      print cov, isolate, path,
            (isolate in mlst_scheme ? mlst_scheme[isolate] : ""),
            (isolate in mlst_st ? mlst_st[isolate] : ""),
            (isolate in mlst_alleles ? mlst_alleles[isolate] : "")
    }
  ' "${MLST_NORM}" "${MANIFEST}" > "${COMBINED_TSV}"

  awk -F'\t' '
    function q(s) { gsub(/"/, "\"\"", s); return "\"" s "\"" }
    BEGIN { OFS="," }
    {
      for (i=1; i<=NF; i++) $i=q($i)
      print
    }
  ' "${COMBINED_TSV}" > "${COMBINED_CSV}"

  awk 'NR==1{next} {print}' "${COMBINED_TSV}" >> "${OVERALL_TSV}"

  echo "[INFO] Wrote:"
  echo "  ${COMBINED_TSV}"
  echo "  ${COMBINED_CSV}"
done

if [[ "${have_any_dir}" -eq 0 ]]; then
  echo "[ERROR] No input directories matched: ${INPUT_GLOB}" >&2
fi

awk -F'\t' '
  function q(s) { gsub(/"/, "\"\"", s); return "\"" s "\"" }
  BEGIN { OFS="," }
  {
    for (i=1; i<=NF; i++) $i=q($i)
    print
  }
' "${OVERALL_TSV}" > "${OVERALL_CSV}"

echo "[INFO] Overall combined outputs:"
echo "  ${OVERALL_TSV}"
echo "  ${OVERALL_CSV}"
                   

  
