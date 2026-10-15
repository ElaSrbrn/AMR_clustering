#!/usr/bin/env python3
"""Add 7-gene MLST sequence types to the cgMLST results.

For each isolate in the cgMLST manifest, takes the mlst ST of the reference assembly
(highest depth >= --min-ref-depth, the same choice as cgmlst_depth.py) and, for comparison,
the ST at every depth. Then labels each pair in cgmlst_pairwise.tsv (reference distances)
with both STs and whether they match.

Inputs:
  --mlst-tables  mlst tables with columns coverage, isolate, mlst_st (07_mlst.sh output)
  --work         cgmlst_work (manifest.tsv)
  --pairwise     cgmlst_pairwise.tsv from cgmlst_depth.py
Outputs:
  <out>_isolates.tsv     id, species, scheme, ref_depth, st_ref, st_by_depth, st_consistent
  <out>_pairs_ref.tsv    scheme, id1, id2, st1, st2, same_st, dist_ref, shared_ref (sorted by distance)

Usage:
  python3 cgmlst_add_st.py --mlst-tables all_coverages_mlst_combined.tsv \
      --work cgmlst_work --pairwise cgmlst_pairwise.tsv --out cgmlst_st
"""
import argparse, csv, re, collections, sys
from pathlib import Path

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--mlst-tables", nargs="+", required=True)
ap.add_argument("--work", default="cgmlst_work")
ap.add_argument("--pairwise", default="cgmlst_pairwise.tsv")
ap.add_argument("--min-ref-depth", type=int, default=60)
ap.add_argument("--out", default="cgmlst_st")
args = ap.parse_args()

num = lambda s: (re.match(r"^isolate_(\d+)", s) or re.match(r"^(\d+)", s) or [None, None])[1]
dnum = lambda c: int(re.sub(r"\D", "", c) or 0)
rd = lambda p: csv.DictReader(open(p, newline=""), delimiter="\t")

# mlst: (id, depth) -> most frequent ST (tables can repeat rows)
votes = collections.defaultdict(collections.Counter)
for t in args.mlst_tables:
    for r in rd(t):
        i = num((r.get("isolate") or "").strip())
        if i:
            votes[(i, dnum(r.get("coverage", "")))][(r.get("mlst_st") or "").strip() or "NA"] += 1
st = {k: c.most_common(1)[0][0] for k, c in votes.items()}

man = list(rd(Path(args.work) / "manifest.tsv"))
iso = {}
for r in man:
    e = iso.setdefault(r["id"], dict(id=r["id"], species=r["species"], scheme=r["scheme"], depths=set()))
    e["depths"].add(int(r["depth"]))

rows = []
for i, e in sorted(iso.items(), key=lambda x: int(x[0])):
    ds = sorted(e["depths"], reverse=True)
    ref = next((d for d in ds if d >= args.min_ref_depth), None)
    by_d = {d: st.get((i, d), "NA") for d in ds}
    numeric = {v for v in by_d.values() if v.isdigit()}
    e["st_ref"] = by_d.get(ref, "NA") if ref else "NA"
    rows.append(dict(id=i, species=e["species"], scheme=e["scheme"], ref_depth=ref or "none",
                     st_ref=e["st_ref"], st_by_depth=";".join(f"{d}x:{by_d[d]}" for d in ds),
                     st_consistent="yes" if len(numeric) == 1 else ("no" if numeric else "NA")))

pairs = {}
if Path(args.pairwise).exists():
    for r in rd(args.pairwise):
        k = (r["scheme"], r["id1"], r["id2"])
        if k in pairs:
            continue
        s1, s2 = iso.get(r["id1"], {}).get("st_ref", "NA"), iso.get(r["id2"], {}).get("st_ref", "NA")
        pairs[k] = dict(scheme=k[0], id1=k[1], id2=k[2], st1=s1, st2=s2,
                        same_st=("yes" if s1 == s2 else "no") if s1.isdigit() and s2.isdigit() else "NA",
                        dist_ref=int(r["dist_ref"]), shared_ref=r["shared_ref"])
else:
    print(f"WARNING: {args.pairwise} not found; writing isolate table only", file=sys.stderr)

def write(name, data):
    if not data:
        return
    with open(name, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(data[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(data)

write(f"{args.out}_isolates.tsv", rows)
write(f"{args.out}_pairs_ref.tsv", sorted(pairs.values(), key=lambda r: (r["scheme"], r["dist_ref"])))

missing = [r["id"] for r in rows if not r["st_ref"].isdigit()]
print(f"isolates: {len(rows)}; numeric ST at reference depth: {len(rows) - len(missing)}", file=sys.stderr)
if missing:
    print(f"no numeric ST at reference depth: {', '.join(missing)}", file=sys.stderr)
print(f"written: {args.out}_isolates.tsv" + (f", {args.out}_pairs_ref.tsv" if pairs else ""), file=sys.stderr)
