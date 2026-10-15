#!/usr/bin/env python3
"""Step 04: summary tables from results/stability_long.tsv.

Outputs: results/summary_by_depth.tsv   counts per assembler x depth x class
         results/stability_rates.tsv    per assembler x depth: n, stable (strict), stable (structural), proportions
         results/stability_matrix.tsv   plasmid x (assembler_depth) -> class   (input for a heatmap)
Definitions:
  stable_strict      class == concordant (includes circularity)
  stable_structural  class in {concordant, concordant_not_circular}
Pairs that were not annotated are excluded from denominators and listed separately.
"""
import sys
from collections import Counter, defaultdict

sys.path.insert(0, "lib")
from common import read_tsv, write_tsv  # noqa: E402

rows = read_tsv("results/stability_long.tsv")
depth_order = {"10x": 0, "20x": 1, "40x": 2, "60x": 3}

cnt = Counter((r["assembler"], r["depth"], r["class"]) for r in rows)
summ = [{"assembler": a, "depth": d, "class": c, "n": n}
        for (a, d, c), n in sorted(cnt.items(), key=lambda x: (x[0][0], depth_order.get(x[0][1], 9), x[0][2]))]
write_tsv("results/summary_by_depth.tsv", summ, ["assembler", "depth", "class", "n"])

groups = defaultdict(list)
for r in rows:
    if r["class"] != "not_annotated":
        groups[(r["assembler"], r["depth"])].append(r["class"])
rates = []
for (a, d), cl in sorted(groups.items(), key=lambda x: (x[0][0], depth_order.get(x[0][1], 9))):
    n = len(cl)
    s1 = sum(c == "concordant" for c in cl)
    s2 = sum(c in ("concordant", "concordant_not_circular") for c in cl)
    rates.append({"assembler": a, "depth": d, "n": n, "stable_strict": s1, "stable_structural": s2,
                  "prop_strict": round(s1 / n, 3), "prop_structural": round(s2 / n, 3)})
write_tsv("results/stability_rates.tsv", rates,
          ["assembler", "depth", "n", "stable_strict", "stable_structural", "prop_strict", "prop_structural"])

cols = sorted({f"{r['assembler']}_{r['depth']}" for r in rows},
              key=lambda x: (x.split("_")[0], depth_order.get(x.split("_")[1], 9)))
mat = defaultdict(dict)
meta = {}
for r in rows:
    mat[r["plasmid_id"]][f"{r['assembler']}_{r['depth']}"] = r["class"]
    meta[r["plasmid_id"]] = {"genes": r["genes"], "ref_length": r["ref_length"], "ref_rep_type": r["ref_rep_type"]}
out = []
for pid in sorted(mat):
    o = {"plasmid_id": pid, **meta[pid]}
    o.update({c: mat[pid].get(c, "") for c in cols})
    out.append(o)
write_tsv("results/stability_matrix.tsv", out, ["plasmid_id", "genes", "ref_length", "ref_rep_type"] + cols)

print("Stability (structural / strict) per assembler and depth:")
for r in rates:
    print(f"  {r['assembler']:10s} {r['depth']:4s}  n={r['n']:3d}  structural={r['prop_structural']:.3f}  strict={r['prop_strict']:.3f}")
na = sum(r["class"] == "not_annotated" for r in rows)
if na:
    print(f"WARNING: {na} pairs not annotated (excluded); see results/stability_long.tsv")
