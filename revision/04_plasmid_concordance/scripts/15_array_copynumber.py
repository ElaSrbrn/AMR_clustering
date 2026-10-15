#!/usr/bin/env python3
"""Step 15 helper: tandem-unit copy number per read, from reads spanning the whole array.
Usage: 15_array_copynumber.py <paf with cg tags> <left_anchor> <right_anchor> <unit_len> <ref_units> <out.tsv>"""
import re, sys
from collections import Counter, defaultdict
paf, LA, RA, UNIT, REF_UNITS, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
def project(ts, te, qs, qe, strand, cigar, pos):
    if not (ts <= pos < te):
        return None
    t, q = ts, 0
    for n, o in CIG.findall(cigar):
        n = int(n)
        if o in "M=X":
            if t + n > pos:
                q += pos - t; break
            t += n; q += n
        elif o == "I":
            q += n
        elif o in "DN":
            if t + n > pos:
                break
            t += n
    return qs + q if strand == "+" else qe - q
alns = defaultdict(list)
for line in open(paf):
    p = line.rstrip("\n").split("\t")
    if len(p) < 12:
        continue
    cg = next((x[5:] for x in p[12:] if x.startswith("cg:Z:")), None)
    if cg is None or any(x == "tp:A:S" for x in p[12:]):
        continue
    alns[p[0]].append((int(p[1]), int(p[2]), int(p[3]), p[4], int(p[7]), int(p[8]), int(p[11]), cg))
rows = []
for read, al in alns.items():
    hits = {}
    for ql, qs, qe, strand, ts, te, mq, cg in al:
        for name, pos in (("L", LA), ("R", RA)):
            q = project(ts, te, qs, qe, strand, cg, pos)
            if q is not None and (name not in hits or mq > hits[name][2]):
                hits[name] = (q, strand, mq)
    if "L" in hits and "R" in hits and hits["L"][1] == hits["R"][1]:
        dq = abs(hits["R"][0] - hits["L"][0]); dr = RA - LA; delta = (dq - dr) / UNIT
        rows.append({"read": read, "read_length": al[0][0], "strand": hits["L"][1], "mapq_left": hits["L"][2],
                     "mapq_right": hits["R"][2], "read_distance": dq, "ref_distance": dr, "delta_units": round(delta, 2),
                     "units_in_read": round(REF_UNITS + delta, 2), "units_rounded": int(round(REF_UNITS + delta))})
cols = ["read", "read_length", "strand", "mapq_left", "mapq_right", "read_distance", "ref_distance", "delta_units", "units_in_read", "units_rounded"]
with open(out, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")
c = Counter(r["units_rounded"] for r in rows)
off = [r for r in rows if abs(r["units_in_read"] - r["units_rounded"]) > 0.15]
print(f"reads spanning both anchors: {len(rows)} (anchor distance in the reference {RA - LA} bp, unit {UNIT} bp, reference units {REF_UNITS})")
for k in sorted(c):
    print(f"  {k} unit(s): {c[k]} reads")
print(f"reads whose estimate is more than 0.15 units from an integer: {len(off)} (indels or other structural differences)")
