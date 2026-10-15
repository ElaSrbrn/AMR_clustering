#!/usr/bin/env python3
"""Step 01c: recompute contigs.tsv (length, circularity) without rerunning AMRFinderPlus/MOB-typer.
Usage: $PY scripts/01c_refresh_contigs.py [assembler]   (default: hybracter)"""
import os, sys
sys.path.insert(0, "lib")
from common import read_tsv, write_tsv, contig_table
asm = sys.argv[1] if len(sys.argv) > 1 else "hybracter"
n = changed = 0
for r in read_tsv("work/manifest.tsv"):
    if r["status"] != "OK" or r["assembler"] != asm:
        continue
    out = f"work/annot/{r['id']}/{asm}_{r['depth']}/contigs.tsv"
    if not os.path.exists(out):
        continue
    old = {c["contig"]: c["circular"] for c in read_tsv(out)}
    rows = contig_table(r["path"])
    write_tsv(out, rows, ["contig", "length", "circular"])
    n += 1
    changed += any(old.get(c["contig"]) != c["circular"] for c in rows)
print(f"refreshed {n} contig tables for {asm}; circularity changed in {changed}")
