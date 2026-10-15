#!/usr/bin/env python3
"""Step 14 helper: summarise PLSDB nearest hits. Usage: 14_plsdb_summary.py <workdir> [plsdb_metadata.tsv]"""
import csv, os, sys
from collections import defaultdict
sys.path.insert(0, "lib")
from common import read_tsv, write_tsv, union_len
d = sys.argv[1]
meta_path = sys.argv[2] if len(sys.argv) > 2 else ""
meta, meta_cols = {}, []
if meta_path and os.path.exists(meta_path):
    with open(meta_path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = rd.fieldnames or []
        acc_col = next((c for c in cols if "ACC" in c.upper()), cols[0] if cols else None)
        meta_cols = [c for c in cols if any(k in c.upper() for k in ("DESCRIPTION", "TAXONOMY_SPECIES", "TAXON_NAME", "SPECIES",
                     "COUNTRY", "LOC_NAME", "COLLECTION_DATE", "CREATE_DATE", "NUCCORE_LENGTH", "LENGTH", "HOST"))][:8]
        for r in rd:
            a = r.get(acc_col, "")
            meta[a] = {c: r.get(c, "") for c in meta_cols}
            meta[a.split(".")[0]] = meta[a]
cov, blocks, qlen = {}, defaultdict(list), {}
paf = os.path.join(d, "aln.paf")
if os.path.exists(paf):
    for line in open(paf):
        p = line.rstrip("\n").split("\t")
        if len(p) < 12 or any(t.startswith("tp:A:") and t != "tp:A:P" for t in p[12:]):
            continue
        q, ql, qs, qe, t, nm, al = p[0], int(p[1]), int(p[2]), int(p[3]), p[5], int(p[9]), int(p[10])
        qlen[q] = ql
        blocks[(q, t)].append((qs, qe, nm, al))
for (q, t), b in blocks.items():
    cov[(q, t)] = (round(100 * union_len([(s, e) for s, e, _, _ in b]) / qlen[q], 1),
                   round(100 * sum(x[2] for x in b) / max(1, sum(x[3] for x in b)), 2))
refs = {r["plasmid_id"]: r for r in read_tsv("results/reference_plasmids.tsv")}
rows = []
for line in open(os.path.join(d, "mash_top.tsv")):
    hit, q, dist, pval, shared = line.rstrip("\n").split("\t")
    r = {"plasmid_id": q, "genes": refs.get(q, {}).get("genes", ""), "ref_length": refs.get(q, {}).get("length", ""),
         "ref_rep_type": refs.get(q, {}).get("rep_type", ""), "plsdb_hit": hit, "mash_distance": dist, "shared_hashes": shared}
    c = cov.get((q, hit))
    r["query_coverage_pct"], r["identity_pct"] = (c if c else ("", ""))
    m = meta.get(hit) or meta.get(hit.split(".")[0]) or {}
    for k in meta_cols:
        r[k] = m.get(k, "")
    rows.append(r)
fields = ["plasmid_id", "genes", "ref_length", "ref_rep_type", "plsdb_hit", "mash_distance", "shared_hashes",
          "query_coverage_pct", "identity_pct"] + meta_cols
write_tsv("results/plsdb_hits.tsv", rows, fields)
seen, nearest = set(), []
for r in rows:
    if r["plasmid_id"] not in seen:
        seen.add(r["plasmid_id"]); nearest.append(r)
write_tsv("results/plsdb_nearest.tsv", nearest, fields)
n = len(nearest)
print(f"written: results/plsdb_hits.tsv ({len(rows)} rows), results/plsdb_nearest.tsv ({n} plasmids)")
print(f"nearest hit with Mash distance <= 0.01: {sum(1 for r in nearest if float(r['mash_distance']) <= 0.01)}/{n}; "
      f"coverage/identity available for {sum(1 for r in nearest if r['query_coverage_pct'] != '')}/{n}")
print("metadata columns used:", ", ".join(meta_cols) or "none")
