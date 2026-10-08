#!/usr/bin/env bash
# run_bart.sh  --  read-based MLST and AMR gene detection with bart v0.1.2 (KMA), one FASTQ at a time.
# Runs on any Linux/macOS machine or cluster node; no job scheduler needed.
#
# Usage:
#   bash run_bart.sh [-t threads] [-m scheme_map.tsv] <reads_dir> [<reads_dir> ...]
#   e.g. bash run_bart.sh -t 8 carblot_pos/10x carblot_pos/20x carblot_neg/10x carblot_neg/20x
#
# Output per sample: <reads_dir>/bart/<sample>/<sample>_mlst.tsv, <sample>_amr.tsv, logs
# Requires: bart with populated MLST and AMR indexes (see README); run inside the bart environment.
# scheme_map.tsv (optional, tab-separated): <sample prefix> <bart scheme name>, to force a scheme.

set -uo pipefail

threads=4; schememap=""
while getopts "t:m:" o; do
    case $o in t) threads=$OPTARG;; m) schememap=$OPTARG;; *) exit 1;; esac
done
shift $((OPTIND-1))
(( $# > 0 )) || { echo "Usage: bash $0 [-t threads] [-m scheme_map.tsv] <reads_dir> ..." >&2; exit 1; }
command -v bart >/dev/null || { echo "ERROR: bart not found (activate the bart environment)" >&2; exit 1; }

n=0; failed=0
shopt -s nullglob
for d in "$@"; do
    [[ -d "$d" ]] || { echo "ERROR: not a directory: $d" >&2; continue; }
    for fq in "$d"/*.fastq "$d"/*.fastq.gz "$d"/*.fq "$d"/*.fq.gz; do
        s="$(basename "$fq")"; s="${s%.gz}"; s="${s%.fastq}"; s="${s%.fq}"
        out="$d/bart/$s"
        mkdir -p "$out"
        # private temp dir: bart writes <sample>.* and choose_scheme.bart into $TMPDIR
        export TMPDIR="$out/tmp_bart"; rm -rf "$TMPDIR"; mkdir -p "$TMPDIR"

        # optional forced MLST scheme (longest matching sample prefix)
        sarg=()
        if [[ -n "$schememap" && -f "$schememap" ]]; then
            sc=$(awk -F'\t' -v s="$s" '$1!="" && (s==$1 || index(s,$1"_")==1) && length($1)>b {b=length($1); x=$2} END{print x}' "$schememap")
            [[ -n "$sc" ]] && sarg=(-s "$sc")
        fi

        # -k is required: without it bart crashes in its cleanup step
        bart "$fq" -r ont -t "$threads" -vv -k "${sarg[@]}" > "$out/${s}_mlst.tsv" 2> "$out/${s}_mlst.log"; mrc=$?
        bart "$fq" -r ont -t "$threads" -amr 90 -k      > "$out/${s}_amr.tsv"  2> "$out/${s}_amr.log";  arc=$?

        rm -rf "$TMPDIR"
        n=$((n+1))
        if (( mrc != 0 || arc != 0 )); then
            failed=$((failed+1)); echo "FAILED $s (mlst=$mrc amr=$arc); see $out/*.log" >&2
        else
            echo "done $s"
        fi
    done
done
echo "$n samples processed, $failed with errors"
