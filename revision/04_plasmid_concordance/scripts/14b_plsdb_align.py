#!/usr/bin/env python3
"""Step 14b: per reference plasmid and its nearest PLSDB hit (Mash): coverage of both plasmids, identity, indels, collinearity.

Input : results/plsdb_nearest.tsv, work/plsdb/nearest_hits.fasta, work/plsdb/nearest.paf
        (minimap2 -c -x asm10 nearest_hits.fasta queries.fasta)
Output: results/plsdb_nearest_alignment.tsv
  query_coverage_pct   share of our plasmid covered by alignment blocks >= 1 kb to the nearest hit
  hit_coverage_pct     share of the PLSDB plasmid covered
  identity_pct         gap-compressed identity, 1 - de (minimap2 de:f tag), weighted by aligned query length;
                       each indel counts once regardless of its size, so large insertions do not lower it
  blast_identity_pct   matches / alignment columns incl. gaps (the previous definition; lowered by large indels)
  max_indel_bp         largest difference between query and target span within one block (proxy for the largest indel)
  n_blocks             alignment blocks >= 1 kb
  collinear            yes if all blocks are on one strand and in the same order (one circular wrap-around allowed)
Diagnostics (stdout): for plasmids without blocks, whether minimap2 reported any line for them and against which target.
"""
import sys
from collections import defaultdict

sys.path.insert(0, "lib")
from common import read_tsv, write_tsv, union_len, read_fasta  # noqa: E402

MIN_BLOCK = 1000


def acc(x):
    return x.split()[0].split(".")[0]


near = {r["plasmid_id"]: r for r in read_tsv("results/plsdb_nearest.tsv")}
hitlen = {acc(n): len(s) for n, h, s in read_fasta("work/plsdb/nearest_hits.fasta")}
blocks, raw = {}, defaultdict(list)
for line in open("work/plsdb/nearest.paf"):
    p = line.rstrip("\n").split("\t")
    if len(p) < 12:
        continue
    q, ql, qs, qe, strand, t, tl, ts, te, nm, al = (p[0], int(p[1]), int(p[2]), int(p[3]), p[4], p[5], int(p[6]),
                                                    int(p[7]), int(p[8]), int(p[9]), int(p[10]))
    tags = dict((x[:2], x[5:]) for x in p[12:] if len(x) > 5 and x[2] == ":" and x[4] == ":")
    raw[q].append((t, qe - qs, tags.get("tp", "")))
    if tags.get("tp", "P") != "P":
        continue
    hit = near.get(q, {}).get("plsdb_hit", "")
    if acc(t) != acc(hit) or qe - qs < MIN_BLOCK:
        continue
    de = float(tags["de"]) if "de" in tags else None
    blocks.setdefault(q, []).append((qs, qe, ts, te, strand, nm, al, ql, tl, de))


def collinear(b):
    if len({x[4] for x in b}) > 1:
        return "no"
    b = sorted(b)
    order = [x[2] for x in b]
    if b[0][4] == "-":
        order = order[::-1]
    breaks = sum(1 for a, c in zip(order, order[1:]) if c < a)   # one wrap-around allowed (circular start offsets)
    return "yes" if breaks <= 1 else "no"


rows = []
for pid, r in near.items():
    b = blocks.get(pid, [])
    out = {"plasmid_id": pid, "genes": r["genes"], "ref_length": r["ref_length"], "ref_rep_type": r["ref_rep_type"],
           "plsdb_hit": r["plsdb_hit"], "hit_length": hitlen.get(acc(r["plsdb_hit"]), ""),
           "mash_distance": r["mash_distance"], "shared_hashes": r["shared_hashes"]}
    if b:
        ql, tl = b[0][7], b[0][8]
        w = [x[1] - x[0] for x in b]
        ide = ("" if any(x[9] is None for x in b)
               else round(100 * (1 - sum(x[9] * wi for x, wi in zip(b, w)) / sum(w)), 2))
        out.update({"query_coverage_pct": round(100 * union_len([(x[0], x[1]) for x in b]) / ql, 1),
                    "hit_coverage_pct": round(100 * union_len([(x[2], x[3]) for x in b]) / tl, 1),
                    "identity_pct": ide,
                    "blast_identity_pct": round(100 * sum(x[5] for x in b) / sum(x[6] for x in b), 2),
                    "max_indel_bp": max(abs((x[3] - x[2]) - (x[1] - x[0])) for x in b),
                    "n_blocks": len(b), "collinear": collinear(b)})
    else:
        out.update({"query_coverage_pct": 0, "hit_coverage_pct": 0, "identity_pct": "", "blast_identity_pct": "",
                    "max_indel_bp": "", "n_blocks": 0, "collinear": ""})
    rows.append(out)

fields = ["plasmid_id", "genes", "ref_length", "ref_rep_type", "plsdb_hit", "hit_length", "mash_distance", "shared_hashes",
          "query_coverage_pct", "hit_coverage_pct", "identity_pct", "blast_identity_pct", "max_indel_bp", "n_blocks", "collinear"]
write_tsv("results/plsdb_nearest_alignment.tsv", rows, fields)


def f(x, k):
    return x[k] != "" and float(x[k])


full = [x for x in rows if x["n_blocks"] and f(x, "query_coverage_pct") >= 95 and f(x, "hit_coverage_pct") >= 95
        and x["identity_pct"] != "" and float(x["identity_pct"]) >= 99.9 and x["collinear"] == "yes"]
print(f"written: results/plsdb_nearest_alignment.tsv ({len(rows)} plasmids)")
print(f">=95% coverage of both plasmids, >=99.9% gap-compressed identity and collinear: {len(full)}/{len(rows)}")
nob = [x for x in rows if x["n_blocks"] == 0]
print(f"no alignment block >= 1 kb to the nearest hit: {len(nob)}")
qnames = set(raw)
print(f"query names in the PAF that are not plasmid_ids in plsdb_nearest.tsv: {sorted(qnames - set(near))[:10]}")
for x in nob:
    lines = raw.get(x["plasmid_id"], [])
    if not lines:
        why = "minimap2 reported no alignment at all for this query"
    else:
        tg = sorted({(acc(t), tp) for t, _, tp in lines})
        longest = max(l for _, l, _ in lines)
        why = f"{len(lines)} PAF lines, targets {tg}, longest {longest} bp"
    print(f"  {x['plasmid_id']}\tnearest {x['plsdb_hit']}\tmash {x['mash_distance']}\t{why}")
