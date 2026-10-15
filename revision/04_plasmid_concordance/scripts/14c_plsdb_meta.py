#!/usr/bin/env python3
"""Step 14c: all PLSDB metadata for the nearest hits (provenance: species, host, country, collection date, BioSample).
Usage: 14c_plsdb_meta.py <plsdb metadata tsv>   ->  results/plsdb_nearest_metadata.tsv (all columns)"""
import csv, sys
sys.path.insert(0, "lib")
from common import read_tsv  # noqa: E402
csv.field_size_limit(10**9)
hits = {r["plsdb_hit"].split(".")[0] for r in read_tsv("results/plsdb_nearest.tsv")}
with open(sys.argv[1]) as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    cols = rd.fieldnames
    acc = next(c for c in cols if "ACC" in c.upper())
    rows = [r for r in rd if r[acc].split(".")[0] in hits]
print("PLSDB metadata columns:", ", ".join(cols))
with open("results/plsdb_nearest_metadata.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, cols, delimiter="\t"); w.writeheader(); w.writerows(rows)
print(f"written: results/plsdb_nearest_metadata.tsv ({len(rows)} of {len(hits)} nearest hits found)")
show = [c for c in cols if any(k in c.upper() for k in ("COUNTRY", "LOC", "GEO", "COLLECTION", "DATE", "HOST", "ISOLATION", "SOURCE", "SPECIES"))]
print("location/date/source columns:", ", ".join(show))
for c in show:
    filled = sum(1 for r in rows if r[c].strip() and r[c].strip().lower() not in ("nan", "na", "missing", "not collected", "unknown"))
    print(f"  {c}: filled for {filled}/{len(rows)}")
