#!/bin/bash
#SBATCH --job-name=pcs_plsdb_aln
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/14b_plsdb_align_%j.out
# Step 14b: align each reference plasmid to its own nearest PLSDB plasmid (pairwise), then summarise.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source ./config.sh
export PATH=${MAMBA_ROOT_PREFIX:-$HOME/mamba}/envs/$ENV_MAP/bin:$PATH
echo "minimap2 $(minimap2 --version)"
[ -s work/plsdb/nearest_hits.fasta ] || { echo "work/plsdb/nearest_hits.fasta missing"; exit 1; }
echo "downloaded hit sequences: $(grep -c '>' work/plsdb/nearest_hits.fasta)"
rm -rf work/plsdb/pairs && mkdir -p work/plsdb/pairs
$PY - <<'EOF'
import sys; sys.path.insert(0, "lib")
from common import read_tsv, read_fasta
near = {r["plasmid_id"]: r["plsdb_hit"] for r in read_tsv("results/plsdb_nearest.tsv")}
q = {n.split()[0]: s for n, h, s in read_fasta("work/plsdb/queries.fasta")}
t = {n.split()[0].split(".")[0]: (n.split()[0], s) for n, h, s in read_fasta("work/plsdb/nearest_hits.fasta")}
missing = []
for pid, hit in near.items():
    a = hit.split(".")[0]
    if pid not in q or a not in t:
        missing.append(pid); continue
    open(f"work/plsdb/pairs/{pid}.q.fa", "w").write(f">{pid}\n{q[pid]}\n")
    n, s = t[a]
    open(f"work/plsdb/pairs/{pid}.t.fa", "w").write(f">{n}\n{s}\n")
print("pairs written:", len(near) - len(missing), "missing:", missing)
EOF
: > work/plsdb/nearest.paf
for qf in work/plsdb/pairs/*.q.fa; do
  p=${qf%.q.fa}
  minimap2 -c -x asm10 -t 2 "$p.t.fa" "$qf" 2>/dev/null >> work/plsdb/nearest.paf
done
echo "PAF lines: $(wc -l < work/plsdb/nearest.paf)"
$PY scripts/14b_plsdb_align.py
