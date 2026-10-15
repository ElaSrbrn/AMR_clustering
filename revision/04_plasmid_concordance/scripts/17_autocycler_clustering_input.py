#!/usr/bin/env python3
"""Step 17a: Mash/pling input for the 60x Autocycler reference plasmids, in the layout that
group_plasmids_mash_then_pling_v2_compat.py expects (same script as for the Flye 40x/60x clustering).
Writes clustering_autocycler60x/annotation_autocycler_60x/{plasmids,summaries}/ from results/reference_plasmids.tsv
and work/compare/<plasmid_id>/ref.fasta. File names carry isolate_<ID>_ for 11_compare_pling_depths.py."""
import os
import sys

sys.path.insert(0, "lib")
from common import read_tsv, read_fasta, write_fasta  # noqa: E402

out = "clustering_autocycler60x/annotation_autocycler_60x"
os.makedirs(f"{out}/plasmids", exist_ok=True)
os.makedirs(f"{out}/summaries", exist_ok=True)
refs = read_tsv("results/reference_plasmids.tsv")
rows, ids = [], {}
for r in refs:
    pid, iso = r["plasmid_id"], r["id"]
    seq = read_fasta(f"work/compare/{pid}/ref.fasta")[0][2]
    fa = os.path.abspath(f"{out}/plasmids/isolate_{iso}_autocycler60x_{pid}.fasta")
    write_fasta(fa, [(f"isolate_{iso}_autocycler60x_{pid}", seq)])
    rep = (r.get("rep_type") or "").strip() or "-"
    rows.append([f"isolate_{iso}", pid, r["genes"], str(len(seq)), rep, "yes", fa])
    ids.setdefault(iso, []).append((pid, rep))
with open(f"{out}/summaries/carbapenemase_summary_autocycler60x.tsv", "w") as fh:
    fh.write("sample\tplasmid_id\tgenes\tlength\treplicon_type\tplasmid_fasta_exported\tplasmid_fasta_path\n")
    for x in rows:
        fh.write("\t".join(x) + "\n")
for iso, lst in ids.items():
    seen = {}
    for pid, rep in lst:
        for g in sorted({t.strip() for t in rep.replace(";", ",").split(",")}):
            if g and g != "-":
                if g in seen:
                    print(f"WARN: isolate {iso} has two {g} plasmids ({seen[g]}, {pid}); depth comparison matches by isolate")
                seen[g] = pid
n = lambda g: sum(1 for x in rows if g in [t.strip() for t in x[4].replace(';', ',').split(',')])
print(f"{len(rows)} reference plasmids written to {out}/ (IncN {n('IncN')}, IncL/M {n('IncL/M')})")
