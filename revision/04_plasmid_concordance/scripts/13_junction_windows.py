#!/usr/bin/env python3
"""Step 13 helper: write one FASTA per junction window for a plasmid (see 13_junction_support.sbatch).
Usage: 13_junction_windows.py <plasmid_id> <outdir>"""
import os
import sys

sys.path.insert(0, "lib")
from common import read_fasta, read_tsv, write_fasta, write_tsv  # noqa: E402

pid, d = sys.argv[1], sys.argv[2]
W = int(os.environ.get("W", 10000))
bps = [b for b in read_tsv("results/breakpoints.tsv") if b["plasmid_id"] == pid and b["assembler"] == "flye"]
if not bps:
    sys.exit(f"{pid}: no Flye breakpoints in results/breakpoints.tsv")
dep = max({b["depth"] for b in bps}, key=lambda x: int(x.rstrip("x")))
bps = [b for b in bps if b["depth"] == dep]
ref = read_fasta(f"work/compare/{pid}/ref.fasta")[0][2]
man = {(r["id"], r["assembler"], r["depth"]): r["path"] for r in read_tsv("work/manifest.tsv")}
contig = bps[0]["primary_contig"]
test = next(s for n, h, s in read_fasta(man[(pid.split("_")[0], "flye", dep)]) if n == contig)
long = {(r["assembler"], r["depth"]): r for r in read_tsv("results/stability_long.tsv") if r["plasmid_id"] == pid}
test_circ = long.get(("flye", dep), {}).get("primary_circular", "false") == "true"


def window(seq, pos, circular):
    """Return (window sequence, 0-based centre in the window)."""
    L = len(seq)
    pos = max(1, min(int(float(pos)), L)) - 1
    w = min(W, L // 2 - 1) if circular else W
    if circular:
        return "".join(seq[i % L] for i in range(pos - w, pos + w)), w
    a, b = max(0, pos - w), min(L, pos + w)
    return seq[a:b], pos - a


wins, rows = {}, []
for i, b in enumerate(bps, 1):
    rk = f"R_{int(float(b['ref_pos']))}"
    tk = f"T_{int(float(b['query_pos']))}" if b.get("query_pos") else ""
    if rk not in wins:
        wins[rk] = window(ref, b["ref_pos"], True)
    if tk and tk not in wins:
        wins[tk] = window(test, b["query_pos"], test_circ)
    rows.append({"plasmid_id": pid, "flye_depth": dep, "breakpoint": i, "class": b["class"], "kind": b["kind"],
                 "side": b["side"], "ref_pos": b["ref_pos"], "query_pos": b.get("query_pos", ""),
                 "ref_window": rk, "test_window": tk})
with open(f"{d}/window_centres.tsv", "w") as fh:
    for k, (s, c) in wins.items():
        write_fasta(f"{d}/win/{k}.fasta", [(k, s)])
        fh.write(f"{k}\t{c}\t{len(s)}\n")
write_tsv(f"{d}/breakpoints.tsv", rows, list(rows[0].keys()))
print(f"  {len(rows)} breakpoints at Flye {dep} ({contig}, circular={test_circ}), {len(wins)} windows; "
      f"reference {len(ref)} bp, Flye contig {len(test)} bp")
