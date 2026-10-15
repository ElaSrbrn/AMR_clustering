#!/usr/bin/env python3
"""Step 12b: summarise step 12 log(s) into results/read_support.tsv (one row per plasmid x sequence).
Low-depth regions (<25% of median, >=500 bp) are split into edge (within EDGE bp of an end: expected for circular
sequences mapped as linear, not informative) and internal (informative).
Usage: $PY scripts/12b_summarise_read_support.py [log ...]   (default: logs/12_read_support_*.out)"""
import glob, os, re, sys
sys.path.insert(0, "lib")
from common import write_tsv
EDGE = int(os.environ.get("EDGE", 5000))
logs = sys.argv[1:] or sorted(glob.glob("logs/12_read_support_*.out"))
rows, cur, reads = [], None, ""
pat = re.compile(r"^\s+(\S+): length (\d+)\s+median depth ([\d.]+)\s+low-depth regions.*?:\s*(.*)$")
for log in logs:
    for line in open(log):
        m = re.match(r"^== (\S+)\s+\(reads: (.+)\)", line)
        if m:
            cur, reads = m.group(1), m.group(2); continue
        m = pat.match(line)
        if not (m and cur):
            continue
        name, L, med, regs = m.group(1), int(m.group(2)), float(m.group(3)), m.group(4).strip()
        edge, internal = [], []
        if regs != "none":
            for r in regs.split(";"):
                mm = re.match(r"\s*(\d+)-(\d+)\s+\(mean ([\d.]+)\)", r)
                if mm:
                    a, b = int(mm.group(1)), int(mm.group(2))
                    (edge if a <= EDGE or b >= L - EDGE else internal).append(f"{a}-{b} (mean {mm.group(3)})")
        rows.append({"plasmid_id": cur, "sequence": "flye_60x" if name.startswith("flye_") else "reference",
                     "contig": name, "length": L, "median_depth": med, "low_internal": "; ".join(internal) or "none",
                     "low_edge": "; ".join(edge) or "none", "n_low_internal": len(internal), "reads": reads,
                     "log": os.path.basename(log)})
seen = {}
for r in rows:
    seen[(r["plasmid_id"], r["sequence"])] = r
rows = sorted(seen.values(), key=lambda r: (r["plasmid_id"], r["sequence"] != "reference"))
for pid in sorted({r["plasmid_id"] for r in rows}):
    ref = next((r for r in rows if r["plasmid_id"] == pid and r["sequence"] == "reference"), None)
    fl = next((r for r in rows if r["plasmid_id"] == pid and r["sequence"] == "flye_60x"), None)
    if ref:
        v = "reference supported" if ref["n_low_internal"] == 0 else "reference has internal low-depth region(s)"
        if fl and fl["length"] > ref["length"] * 1.001:
            v += ("; Flye contig longer than the reference and fully read-supported (possible alternative configuration)"
                  if fl["n_low_internal"] == 0 else "; Flye contig longer than the reference, partly unsupported")
        elif fl and fl["length"] < ref["length"] * 0.999:
            v += "; Flye contig shorter (reference sequence missing from Flye)"
        for r in (ref, fl):
            if r: r["interpretation"] = v
os.makedirs("results", exist_ok=True)
write_tsv("results/read_support.tsv", rows, ["plasmid_id", "sequence", "contig", "length", "median_depth",
          "n_low_internal", "low_internal", "low_edge", "interpretation", "reads", "log"])
for r in rows:
    if r["sequence"] == "reference":
        print(f"{r['plasmid_id']:9s} {r['interpretation']}")
print(f"written: results/read_support.tsv ({len(rows)} rows)")
