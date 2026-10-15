#!/usr/bin/env python3
"""Step 05: propose plasmids for in-depth analysis and for synteny figures.

The script only PROPOSES. It writes results/proposal_indepth.tsv and results/proposal_synteny.tsv
with the reason for each plasmid and an 'include' column (yes / review).
Copy them to selection/indepth.tsv and selection/synteny.tsv, edit 'include', and steps 06-08
read only rows with include == yes from the selection/ copies (so re-running 05 never overwrites your edits).

In-depth (ISfinder, Bakta, repeats, breakpoint mapping), proposed if any of:
  A  not structurally stable at 40x or 60x with Flye or Hybracter
  B  Flye and Hybracter disagree (one stable, one not) at the same depth >= 20x
  C  more than one gene copy in the reference (possible duplication)
Synteny (collaborator, pyGenomeViz), proposed:
  1  one representative per failure class (highest depth at which the class occurs, Flye preferred)
  2  small plasmids (< SMALL_PLASMID bp)
  3  plasmids listed in manual_include_synteny.tsv (e.g. Mash/DCJ cluster medoids)
"""
import os
import sys
from collections import defaultdict

sys.path.insert(0, "lib")
from common import read_tsv, write_tsv  # noqa: E402

SMALL = int(os.environ.get("SMALL_PLASMID", 20000))
STABLE = {"concordant", "concordant_not_circular"}
SKIP = {"not_annotated"}
depth_rank = {"10x": 0, "20x": 1, "40x": 2, "60x": 3}

rows = [r for r in read_tsv("results/stability_long.tsv") if r["class"] not in SKIP]
refs = {r["plasmid_id"]: r for r in read_tsv("results/reference_plasmids.tsv")}
by_p = defaultdict(dict)
for r in rows:
    by_p[r["plasmid_id"]][(r["assembler"], r["depth"])] = r["class"]

ind = []
for pid, d in sorted(by_p.items()):
    reasons = []
    bad = [f"{a}_{dp}:{c}" for (a, dp), c in d.items() if dp in ("40x", "60x") and c not in STABLE]
    if bad:
        reasons.append("A " + ";".join(sorted(bad)))
    for dp in ("20x", "40x", "60x"):
        f, h = d.get(("flye", dp)), d.get(("hybracter", dp))
        if f and h and ((f in STABLE) != (h in STABLE)):
            reasons.append(f"B {dp} flye={f} hybracter={h}")
    if int(refs.get(pid, {}).get("n_gene_copies", 1) or 1) > 1:
        reasons.append("C multiple gene copies in reference")
    if reasons:
        ind.append({"plasmid_id": pid, "genes": refs[pid]["genes"], "length": refs[pid]["length"],
                    "rep_type": refs[pid]["rep_type"], "reason": " | ".join(reasons), "include": "yes"})

syn = {}
best = {}
for r in rows:
    c = r["class"]
    if c in STABLE or c == "gene_absent":
        continue
    key = (depth_rank.get(r["depth"], -1), r["assembler"] == "flye")
    if c not in best or key > best[c][0]:
        best[c] = (key, r)
for c, (_, r) in best.items():
    syn.setdefault(r["plasmid_id"], []).append(f"1 representative of '{c}' ({r['assembler']} {r['depth']})")
for pid, ref in refs.items():
    if int(ref["length"]) < SMALL:
        syn.setdefault(pid, []).append(f"2 small plasmid ({ref['length']} bp)")
for m in read_tsv("manual_include_synteny.tsv"):
    if m.get("plasmid_id"):
        syn.setdefault(m["plasmid_id"], []).append("3 manual: " + m.get("reason", ""))
syn_rows = [{"plasmid_id": pid, "genes": refs.get(pid, {}).get("genes", ""), "length": refs.get(pid, {}).get("length", ""),
             "rep_type": refs.get(pid, {}).get("rep_type", ""), "reason": " | ".join(rs), "include": "yes"}
            for pid, rs in sorted(syn.items())]

f = ["plasmid_id", "genes", "length", "rep_type", "reason", "include"]
write_tsv("results/proposal_indepth.tsv", ind, f)
write_tsv("results/proposal_synteny.tsv", syn_rows, f)
os.makedirs("selection", exist_ok=True)
print(f"Proposed for in-depth analysis: {len(ind)}  -> results/proposal_indepth.tsv")
print(f"Proposed for synteny figures:  {len(syn_rows)}  -> results/proposal_synteny.tsv")
print("Review, then copy to selection/indepth.tsv and selection/synteny.tsv and edit the 'include' column.")
