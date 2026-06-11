#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=64G
#SBATCH -t 12:00:00

# ====== CONFIG ======
G="5m"
SEED=24
OUTBASE="../"
TMPROOT="${OUTBASE}/tmp_rasusa"
COVS=(100 80 60 40 30 20 10)
# ====================

mkdir -p "$TMPROOT"
for c in "${COVS[@]}"; do
  mkdir -p "${OUTBASE}/${c}x"
done

REPORT_ERR="${OUTBASE}/rasusa_errors.txt"
REPORT_LOW="${OUTBASE}/low_coverage.txt"
REPORT_RUN="${OUTBASE}/rasusa_runlog.txt"
REPORT_QC="${OUTBASE}/rasusa_qc.tsv"
REPORT_INC="${OUTBASE}/included_nested.tsv"
REPORT_SUMMARY="${OUTBASE}/coverage_nested_summary.txt"

: > "$REPORT_ERR"
: > "$REPORT_LOW"
: > "$REPORT_RUN"
printf "sample\tstage\treads\tbases\tcov_x\tfile\n" > "$REPORT_QC"
printf "sample\toriginal_file\tinput_actual_cov\tmax_nested_cov\toutput_file\n" > "$REPORT_INC"

is_gzip() {
  head -c2 -- "$1" 2>/dev/null | od -An -t x1 | tr -d ' \n' | grep -q '^1f8b$'
}

normalize_name() {
  local bn="$1" cov="$2" stem="$bn"
  stem="${stem%.fastq.gz}"
  stem="${stem%.fastq}"
  [[ "$stem" == *"_filtered" ]] && stem="${stem%_filtered}"
  [[ "$stem" == *"_merged" ]] && stem="${stem%_merged}"
  printf "%s_%s.fastq" "$stem" "$cov"
}

parse_gbp() {
  awk -v s="$1" '
    BEGIN{
      unit=""
      if (s ~ /[kK]$/) unit="k"
      else if (s ~ /[mM]$/) unit="m"
      else if (s ~ /[gG]$/) unit="g"
      g=s
      sub(/[kKmMgG]$/, "", g)
      val=g+0.0
      mult=(unit=="k"?1e3:(unit=="m"?1e6:(unit=="g"?1e9:1)))
      printf("%.0f", val*mult)
    }'
}

cov_fmt() {
  awk -v b="$1" -v g="$2" 'BEGIN{ if(g>0) printf("%.2f", b/g); else print "0.00" }'
}

total_bases() {
  local f="$1"
  if is_gzip "$f"; then
    gzip -cd -- "$f" 2>/dev/null | awk 'NR%4==2{n+=length($0)} END{print n+0}'
  else
    awk 'NR%4==2{n+=length($0)} END{print n+0}' "$f"
  fi
}

total_reads() {
  local f="$1"
  if is_gzip "$f"; then
    gzip -cd -- "$f" 2>/dev/null | awk 'END{print NR/4}'
  else
    awk 'END{print NR/4}' "$f"
  fi
}

prepare_src_fastq() {
  local f="$1" tmpdir="$2" out
  if is_gzip "$f"; then
    out="${tmpdir}/$(basename "${f%.gz}").plain.fastq"
    gzip -cd -- "$f" > "$out" 2>>"$REPORT_ERR"
    printf "%s\n" "$out"
  else
    printf "%s\n" "$f"
  fi
}

gen_cov() {
  local cov="$1" src="$2" out="$3"
  rasusa reads -g "$G" -c "$cov" -s "$SEED" "$src" > "$out"
}

qc_fastq() {
  local label="$1" stage="$2" f="$3"
  local reads bases cov
  reads=$(awk 'END{print NR/4}' "$f")
  bases=$(awk 'NR%4==2{n+=length($0)} END{print n+0}' "$f")
  cov=$(cov_fmt "$bases" "$G_BP")
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$label" "$stage" "$reads" "$bases" "$cov" "$f" >> "$REPORT_QC"
}

G_BP=$(parse_gbp "$G")

n_total=0
n_any=0

shopt -s nullglob

for f in *.fastq *.fastq.gz; do
  ((n_total++))

  bn="$(basename "$f")"
  label="${bn%.fastq.gz}"
  label="${label%.fastq}"

  echo "==> $bn" | tee -a "$REPORT_RUN"

  bases=$(total_bases "$f") || {
    echo "ERROR: $bn : failed to read bases" | tee -a "$REPORT_ERR"
    continue
  }

  reads=$(total_reads "$f") || {
    echo "ERROR: $bn : failed to read reads" | tee -a "$REPORT_ERR"
    continue
  }

  actual_cov="$(cov_fmt "$bases" "$G_BP")"
  echo "input_reads=${reads} ; input_bases=${bases} ; actual_cov=${actual_cov}x" | tee -a "$REPORT_RUN"

  tmpdir="$(mktemp -d -p "$TMPROOT" rasusa.XXXXXX)"
  src=$(prepare_src_fastq "$f" "$tmpdir") || {
    echo "ERROR: $bn : extraction failed" | tee -a "$REPORT_ERR"
    rm -rf "$tmpdir"
    continue
  }

  prev="$src"
  made_any=0
  max_cov="none"

  for cov in "${COVS[@]}"; do
    if (( bases < cov * G_BP )); then
      echo -e "${bn}\tactual_cov=${actual_cov}x\tskipped=${cov}x\treason=total_cov_<${cov}x" >> "$REPORT_LOW"
      continue
    fi

    out="${OUTBASE}/${cov}x/$(normalize_name "$bn" "${cov}x")"

    if ! gen_cov "$cov" "$prev" "$out"; then
      echo "ERROR: $bn : rasusa failure during ${cov}x generation" | tee -a "$REPORT_ERR"
      echo -e "${bn}\tactual_cov=${actual_cov}x\tskipped=${cov}x\treason=rasusa_failure" >> "$REPORT_LOW"
      break
    fi

    qc_fastq "$label" "${cov}x" "$out"
    printf "%s\t%s\t%sx\t%sx\t%s\n" "$label" "$bn" "$actual_cov" "$cov" "$out" >> "$REPORT_INC"

    prev="$out"
    made_any=1
    [[ "$max_cov" == "none" ]] && max_cov="${cov}x"
  done

  (( made_any )) && ((n_any++))

  echo "OK: $bn done" | tee -a "$REPORT_RUN"
  rm -rf "$tmpdir"
done

{
  echo "Nested coverage summary"
  echo "======================="
  echo "Total input isolates: $n_total"
  echo "Generated at least one nested subset: $n_any"
  echo
  echo "Coverage levels:"
  printf " - %sx\n" "${COVS[@]}"
  echo
  echo "Files:"
  echo " - Included: $REPORT_INC"
  echo " - QC table: $REPORT_QC"
  echo " - Errors:   $REPORT_ERR"
  echo " - Skipped:  $REPORT_LOW"
} > "$REPORT_SUMMARY"

echo "Reports:"
echo " - Errors:   $REPORT_ERR"
echo " - Skipped:  $REPORT_LOW"
echo " - Run log:  $REPORT_RUN"
echo " - QC table: $REPORT_QC"
echo " - Included: $REPORT_INC"
echo " - Summary:  $REPORT_SUMMARY"

