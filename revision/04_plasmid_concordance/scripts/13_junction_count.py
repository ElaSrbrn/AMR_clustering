#!/usr/bin/env python3
"""Step 13 helper: count junction-spanning reads per window and classify each breakpoint.
Usage: 13_junction_count.py <outdir>   (prints rows for results/junction_support.tsv)"""
import csv
import os
import re
import sys

d = sys.argv[1]
OV = int(os.environ.get("OVERHANG", 3000))
MINID = float(os.environ.get("MINID", 0.80))
BIG = int(os.environ.get("BIGINDEL", 100))
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
centre = {}
for line in open(f"{d}/window_centres.tsv"):
    k, c, _ = line.rstrip("\n").split("\t")
    centre[k] = int(c)


def big_indel_near(ts, cigar, c):
    """True if the alignment has an insertion/deletion >= BIG within OV of the window centre c (target coordinates)."""
    t = ts
    for n, o in CIG.findall(cigar):
        n = int(n)
        if o in "M=X":
            t += n
        elif o == "I":
            if n >= BIG and c - OV <= t <= c + OV:
                return True
        elif o in "DN":
            if n >= BIG and t <= c + OV and t + n >= c - OV:
                return True
            t += n
    return False


spanning = {}
for k, c in centre.items():
    reads = set()
    path = f"{d}/paf/{k}.paf"
    if os.path.exists(path):
        for line in open(path):
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            ts, te, nm, al = int(p[7]), int(p[8]), int(p[9]), int(p[10])
            if not (ts <= c - OV and te >= c + OV) or nm / max(al, 1) < MINID:
                continue
            cg = next((x[5:] for x in p[12:] if x.startswith("cg:Z:")), "")
            if cg and big_indel_near(ts, cg, c):
                continue
            reads.add(p[0])
    spanning[k] = reads

for r in csv.DictReader(open(f"{d}/breakpoints.tsv"), delimiter="\t"):
    R = spanning.get(r["ref_window"], set())
    T = spanning.get(r["test_window"], set()) if r["test_window"] else set()
    both = R & T
    nr, nt, nb = len(R - T), len(T - R), len(both)
    tot = nr + nt
    frac = nt / tot if tot else None
    if not r["test_window"]:
        verdict = "no Flye position for this breakpoint"
    elif tot < 5:
        verdict = "too few informative reads"
    elif frac < 0.05:
        verdict = "Flye junction unsupported (assembly error)"
    elif frac > 0.95:
        verdict = "reference junction unsupported (reference error)"
    else:
        verdict = "both configurations supported (heterogeneity)"
    if nb > tot:
        verdict += "; windows largely identical (repeat-mediated junction), interpret with care"
    print("\t".join([r["plasmid_id"], r["flye_depth"], r["breakpoint"], r["class"], r["kind"], r["side"], r["ref_pos"],
                     r["query_pos"], str(nr), str(nt), str(nb), "NA" if frac is None else f"{frac:.3f}", verdict]))
