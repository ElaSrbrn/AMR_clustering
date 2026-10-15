#!/usr/bin/env python3
"""Step 06a: prepare in-depth analysis folders for the selected plasmids.

Reads selection/indepth.tsv and selection/synteny.tsv (rows with include == yes; the union of both is used,
because the synteny packages also need the Bakta and IS annotation).
For each plasmid writes work/indepth/<plasmid_id>/:
  ref.fasta            reference plasmid
  test_contigs.fasta   primary contigs of all NON-stable comparisons, named <assembler>_<depth>__<contig>
and the list work/indepth_list.txt used as the array index for 06_indepth.sbatch.
"""
import os
import sys

sys.path.insert(0, "lib")
from common import read_tsv, read_fasta, write_fasta  # noqa: E402

STABLE = {"concordant", "concordant_not_circular"}
sel = []
for f in ("selection/indepth.tsv", "selection/synteny.tsv"):
    if not os.path.exists(f):
        print(f"note: {f} not found (copy from results/proposal_*.tsv and edit 'include')")
        continue
    sel += [r["plasmid_id"] for r in read_tsv(f) if r.get("include", "").strip().lower() == "yes"]
sel = sorted(set(sel))
pairs = {(p["plasmid_id"], p["assembler"], p["depth"]): p for p in read_tsv("work/pairs.tsv")}
long = read_tsv("results/stability_long.tsv")
cache = {}
for pid in sel:
    d = f"work/indepth/{pid}"
    os.makedirs(d, exist_ok=True)
    ref = read_fasta(f"work/compare/{pid}/ref.fasta")
    write_fasta(f"{d}/ref.fasta", [(pid, ref[0][2])])
    tests = []
    for r in long:
        if r["plasmid_id"] != pid or r["class"] in STABLE or not r.get("primary_contig"):
            continue
        fa = pairs[(pid, r["assembler"], r["depth"])]["test_fasta"]
        if fa not in cache:
            cache[fa] = {n: s for n, _, s in read_fasta(fa)}
        seq = cache[fa].get(r["primary_contig"])
        if seq and len(seq) < 1000000:
            tests.append((f"{r['assembler']}_{r['depth']}__{r['primary_contig']}", seq))
    write_fasta(f"{d}/test_contigs.fasta", tests)
with open("work/indepth_list.txt", "w") as f:
    f.write("\n".join(sel) + ("\n" if sel else ""))
print(f"{len(sel)} plasmids prepared in work/indepth/; array size for 06_indepth.sbatch: {len(sel)}")
