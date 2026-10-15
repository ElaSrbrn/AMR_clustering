#!/bin/bash
#SBATCH --job-name=pcs_depth
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --exclude=cpusrv08
#SBATCH --output=logs/10_actual_depth_%j.out
#
# Step 10 (editor: depth definition with a fixed 5 Mb genome size).
# Actual depth = total bases of each read subset / length of that isolate's reference assembly
# (Autocycler 60x consensus, chromosome + plasmids). Isolates without a 60x reference use the
# Flye 40x assembly length as the denominator (flagged in the 'denominator' column).
# Output: results/actual_depth.tsv  (one row per isolate x nominal depth)
#         results/actual_depth_summary.tsv  (per nominal depth and species)
# Usage (from the pipeline folder): sbatch scripts/10_actual_depth.sbatch
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source config.sh
make_jobtmp
eval "$(micromamba shell hook --shell bash)"
micromamba activate "$ENV_MAP"
command -v seqkit >/dev/null || export PATH="${MAMBA_ROOT_PREFIX:-$HOME/mamba}/envs/qc_env/bin:$PATH"
command -v seqkit >/dev/null || { echo "seqkit not found in $ENV_MAP or qc_env"; exit 1; }

OUT=work/actual_depth; mkdir -p $OUT results
for D in $DEPTHS; do
  seqkit stats -T -j 8 $BASE/reads_${D}_S1/*.fastq* > $OUT/reads_${D}.tsv
done
# reference lengths: Autocycler (all S1 isolates with a 60x folder), Flye 40x as fallback
: > $OUT/ref_lengths.tsv
while read -r d; do
  id=$(basename "$d" | cut -d_ -f2)
  f="$d/autocycler_out/consensus_assembly.fasta"
  [ -s "$f" ] && printf '%s\tautocycler_60x\t%s\n' "$id" "$(seqkit stats -T "$f" | awk 'NR==2{print $5}')" >> $OUT/ref_lengths.tsv
done < "$ISOLATE_DIRS"
for id in 5371 5668 6804 7418; do
  g=${FLYE_GLOB//\{ID\}/$id}; g=${g//\{DEPTH\}/40x}
  f=$(ls $g 2>/dev/null | head -1)
  [ -n "$f" ] && printf '%s\tflye_40x\t%s\n' "$id" "$(seqkit stats -T "$f" | awk 'NR==2{print $5}')" >> $OUT/ref_lengths.tsv
done

$PY - <<'PY'
import csv, glob, os, re, statistics, sys
sys.path.insert(0, "lib")
from common import write_tsv
ref = {}
for l in open("work/actual_depth/ref_lengths.tsv"):
    i, src, L = l.rstrip("\n").split("\t"); ref[i] = (src, int(L))
species = {}
if os.path.exists("s1_species.tsv"):
    for l in open("s1_species.tsv"):
        i, s = l.rstrip("\n").split("\t")[:2]; species[i] = s
rows = []
for f in sorted(glob.glob("work/actual_depth/reads_*.tsv")):
    nom = re.search(r"reads_(\d+x)\.tsv", f).group(1)
    for r in csv.DictReader(open(f), delimiter="\t"):
        m = re.search(r"isolate_(\d{4})_", os.path.basename(r["file"]))
        if not m: continue
        i = m.group(1); bases = int(r["sum_len"])
        src, L = ref.get(i, ("NA", 0))
        rows.append({"id": i, "species": species.get(i, ""), "nominal_depth": nom, "bases": bases,
                     "denominator": src, "reference_length": L,
                     "actual_depth": round(bases / L, 1) if L else "",
                     "ratio_actual_nominal": round(bases / L / int(nom[:-1]), 3) if L else ""})
write_tsv("results/actual_depth.tsv", rows, ["id", "species", "nominal_depth", "bases", "denominator",
          "reference_length", "actual_depth", "ratio_actual_nominal"])
summ = []
groups = {}
for r in rows:
    if r["actual_depth"] == "": continue
    for key in ((r["nominal_depth"], "all"), (r["nominal_depth"], r["species"] or "unknown")):
        groups.setdefault(key, []).append(r["actual_depth"])
for (nom, sp), v in sorted(groups.items(), key=lambda x: (int(x[0][0][:-1]), x[0][1] != "all", x[0][1])):
    v = sorted(v)
    q = statistics.quantiles(v, n=4) if len(v) > 1 else [v[0]] * 3
    summ.append({"nominal_depth": nom, "species": sp, "n": len(v), "median": round(statistics.median(v), 1),
                 "q1": round(q[0], 1), "q3": round(q[2], 1), "min": v[0], "max": v[-1]})
write_tsv("results/actual_depth_summary.tsv", summ, ["nominal_depth", "species", "n", "median", "q1", "q3", "min", "max"])
for s in summ:
    if s["species"] == "all":
        print(f"{s['nominal_depth']:>4} nominal: actual median {s['median']} (IQR {s['q1']}-{s['q3']}, range {s['min']}-{s['max']}), n={s['n']}")
PY
