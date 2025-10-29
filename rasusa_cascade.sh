#!/usr/bin/env bash
set -euo pipefail

# ====== CONFIG ======
G="5m"          # genome size; supports 5m, 4.8m, 5000000, etc.
SEED=42
OUTBASE="../"   # where 100x/70x/40x/10x live
TMPROOT="${OUTBASE}/tmp_rasusa"  # avoid /tmp filling up
# ====================

mkdir -p "${OUTBASE}/100x" "${OUTBASE}/70x" "${OUTBASE}/40x" "${OUTBASE}/10x" "$TMPROOT"

REPORT_ERR="${OUTBASE}/rasusa_errors.txt"
REPORT_LOW="${OUTBASE}/low_coverage.txt"
REPORT_RUN="${OUTBASE}/rasusa_runlog.txt"
: > "$REPORT_ERR"; : > "$REPORT_LOW"; : > "$REPORT_RUN"

# --- helpers ---
is_gzip(){ head -c2 -- "$1" 2>/dev/null | od -An -t x1 | tr -d ' \n' | grep -q '^1f8b$'; }

normalize_name(){
  local bn="$1" cov="$2" stem="$bn"
  stem="${stem%.fastq.gz}"; stem="${stem%.fastq}"
  [[ "$stem" == *"_filtered" ]] && stem="${stem%_filtered}"
  [[ "$stem" == *"_merged"   ]] && stem="${stem%_merged}"
  printf "%s_%s.fastq" "$stem" "$cov"
}

# Parse G to integer bp (handles decimals + k/m/g)
parse_gbp(){
  awk -v s="$1" '
    BEGIN{
      unit=""; if (s ~ /[kK]$/) unit="k"; else if (s ~ /[mM]$/) unit="m"; else if (s ~ /[gG]$/) unit="g";
      g=s; sub(/[kKmMgG]$/,"",g); val=g+0.0;
      mult=(unit=="k"?1e3:(unit=="m"?1e6:(unit=="g"?1e9:1)));
      printf("%.0f", val*mult);
    }'
}
G_BP=$(parse_gbp "$G")

cov_fmt(){ awk -v b="$1" -v g="$2" 'BEGIN{ if(g>0) printf("%.2f", b/g); else print "0.00" }'; }

total_bases(){
  local f="$1"
  if is_gzip "$f"; then
    gzip -cd -- "$f" 2>/dev/null | awk 'NR%4==2{n+=length($0)} END{print n+0}'
  else
    awk 'NR%4==2{n+=length($0)} END{print n+0}' "$f"
  fi
}

prepare_src_fastq(){
  # Ensure plain FASTQ on disk for rasusa (decompress if gz or gz-with-trailer)
  local f="$1" tmpdir="$2" out
  if is_gzip "$f"; then
    out="${tmpdir}/$(basename "${f%.*}").plain.fastq"
    if ! gzip -cd -- "$f" > "$out" 2>>"$REPORT_ERR"; then
      return 1
    fi
    printf "%s\n" "$out"
  else
    printf "%s\n" "$f"
  fi
}

gen_cov(){ local cov="$1" src="$2" out="$3"; rasusa reads -g "$G" -c "$cov" -s "$SEED" "$src" > "$out"; }

# --- main ---
shopt -s nullglob
for f in *.fastq *.fastq.gz; do
  [[ -e "$f" ]] || continue
  bn="$(basename "$f")"
  echo "==> $bn" | tee -a "$REPORT_RUN"

  # count bases (FIX: no stray '--' before filename)
  if ! bases=$(total_bases "$f" 2>/dev/null); then
    echo "ERROR: $bn : failed to read bases (gzip/corruption?)" | tee -a "$REPORT_ERR"
    echo -e "${bn}\tactual_cov=NA\tgenerated=none\tskipped=100x 70x 40x 10x\treason=parse_error" >> "$REPORT_LOW"
    continue
  fi
  actual_cov="$(cov_fmt "$bases" "$G_BP")"
  echo "bases=${bases} ; actual_cov=${actual_cov}x" >> "$REPORT_RUN"

  # feasibility
  can100=$(( bases >= 100 * G_BP ? 1 : 0 ))
  can70=$((  bases >=  70 * G_BP ? 1 : 0 ))
  can40=$((  bases >=  40 * G_BP ? 1 : 0 ))
  can10=$((  bases >=  10 * G_BP ? 1 : 0 ))

  # tmp dir (FIX: avoid /tmp)
  tmpdir=$(mktemp -d -p "$TMPROOT")
  # shellcheck disable=SC2064
  trap "rm -rf '$tmpdir'" RETURN

  if ! src=$(prepare_src_fastq "$f" "$tmpdir"); then
    echo "ERROR: $bn : gzip extraction failed" | tee -a "$REPORT_ERR"
    echo -e "${bn}\tactual_cov=${actual_cov}x\tgenerated=none\tskipped=100x 70x 40x 10x\treason=extract_failed" >> "$REPORT_LOW"
    continue
  fi

  out100="${OUTBASE}/100x/$(normalize_name "$bn" 100x)"
  out70="${OUTBASE}/70x/$(normalize_name "$bn" 70x)"
  out40="${OUTBASE}/40x/$(normalize_name "$bn" 40x)"
  out10="${OUTBASE}/10x/$(normalize_name "$bn" 10x)"

  generated=(); skipped=()

  if (( can100 )); then
    { gen_cov 100 "$src" "$out100" \
      && gen_cov 70  "$out100" "$out70" \
      && gen_cov 40  "$out70"  "$out40" \
      && gen_cov 10  "$out40"  "$out10"; } \
    || { echo "ERROR: $bn : rasusa failure in 100/70/40/10 cascade" | tee -a "$REPORT_ERR"; echo -e "${bn}\tactual_cov=${actual_cov}x\tgenerated=${generated[*]:-none}\tskipped=unknown\treason=rasusa_failure" >> "$REPORT_LOW"; continue; }
    generated+=(100x 70x 40x 10x)

  elif (( can70 )); then
    { gen_cov 70 "$src" "$out70" \
      && gen_cov 40 "$out70" "$out40" \
      && gen_cov 10 "$out40" "$out10"; } \
    || { echo "ERROR: $bn : rasusa failure in 70/40/10 cascade" | tee -a "$REPORT_ERR"; echo -e "${bn}\tactual_cov=${actual_cov}x\tgenerated=${generated[*]:-none}\tskipped=unknown\treason=rasusa_failure" >> "$REPORT_LOW"; continue; }
    generated+=(70x 40x 10x); skipped+=(100x)

  elif (( can40 )); then
    { gen_cov 40 "$src" "$out40" \
      && gen_cov 10 "$out40" "$out10"; } \
    || { echo "ERROR: $bn : rasusa failure in 40/10 cascade" | tee -a "$REPORT_ERR"; echo -e "${bn}\tactual_cov=${actual_cov}x\tgenerated=${generated[*]:-none}\tskipped=unknown\treason=rasusa_failure" >> "$REPORT_LOW"; continue; }
    generated+=(40x 10x); skipped+=(100x 70x)

  elif (( can10 )); then
    { gen_cov 10 "$src" "$out10"; } \
    || { echo "ERROR: $bn : rasusa failure at 10x" | tee -a "$REPORT_ERR"; echo -e "${bn}\tactual_cov=${actual_cov}x\tgenerated=${generated[*]:-none}\tskipped=unknown\treason=rasusa_failure" >> "$REPORT_LOW"; continue; }
    generated+=(10x); skipped+=(100x 70x 40x)

  else
    skipped+=(100x 70x 40x 10x)
  fi

  # log skipped tiers w/ actual coverage
  if (( ${#skipped[@]} > 0 )); then
    reason=$([ $can10 -eq 0 ] && echo "total_cov_<10x" || echo "not_enough_for_some")
    echo -e "${bn}\tactual_cov=${actual_cov}x\tgenerated=${generated[*]:-none}\tskipped=${skipped[*]}\treason=${reason}" >> "$REPORT_LOW"
  fi

  echo "OK: $bn done" >> "$REPORT_RUN"
  rm -rf "$tmpdir"
done

echo "Reports:"
echo " - Errors:        $REPORT_ERR"
echo " - Low coverage:  $REPORT_LOW"
echo " - Run log:       $REPORT_RUN"
