#!/usr/bin/env python3
import csv, glob, itertools, os, re, sys
A, B = sys.argv[1], sys.argv[2]
groups = sys.argv[3:] or ["IncN", "IncL_M"]
ISO = re.compile(r"isolate_(\d{4})")
def iso(name):
    m = ISO.search(os.path.basename(name)); return m.group(1) if m else os.path.basename(name)
def load(annot, grp):
    gdir = os.path.join(annot, "mash_then_pling", grp)
    kept = {iso(p) for p in glob.glob(os.path.join(gdir, "mash_within_*", "*"))}
    allp = {iso(p) for p in glob.glob(os.path.join(gdir, "all_plasmids", "*"))}
    typing = {}
    cands = [os.path.join(gdir, "pling", "typing.tsv")] + sorted(glob.glob(os.path.join(gdir, "pling", "**", "typing.tsv"), recursive=True))
    for t in [c for c in cands if os.path.exists(c)]:
        with open(t) as fh:
            rd = csv.reader(fh, delimiter="\t"); next(rd)
            for row in rd:
                if len(row) >= 2: typing[iso(row[0])] = row[1]
        break
    return allp, kept, typing
def comm(label):
    m = re.match(r"(community_\d+)", label); return m.group(1) if m else label
rows = []
for g in groups:
    a_all, a_kept, a_typ = load(A, g); b_all, b_kept, b_typ = load(B, g)
    shared = sorted(set(a_typ) & set(b_typ)); pairs = list(itertools.combinations(shared, 2))
    ssa = {p for p in pairs if a_typ[p[0]] == a_typ[p[1]]}; ssb = {p for p in pairs if b_typ[p[0]] == b_typ[p[1]]}
    sca = {p for p in pairs if comm(a_typ[p[0]]) == comm(a_typ[p[1]])}; scb = {p for p in pairs if comm(b_typ[p[0]]) == comm(b_typ[p[1]])}
    print(f"== {g}")
    print(f"  plasmids: {A} {len(a_all)}, {B} {len(b_all)}")
    print(f"  kept by Mash: {A} {len(a_kept)}, {B} {len(b_kept)}; both {len(a_kept & b_kept)}; only {A}: {sorted(a_kept - b_kept)}; only {B}: {sorted(b_kept - a_kept)}")
    print(f"  pling-typed at both depths: {len(shared)}")
    print(f"  communities: {A} {len({comm(v) for v in a_typ.values()})}, {B} {len({comm(v) for v in b_typ.values()})}; subcommunities: {A} {len(set(a_typ.values()))}, {B} {len(set(b_typ.values()))}")
    if pairs:
        ag_s = sum((p in ssa) == (p in ssb) for p in pairs); ag_c = sum((p in sca) == (p in scb) for p in pairs)
        print(f"  pairwise agreement, subcommunity: {ag_s}/{len(pairs)} ({100*ag_s/len(pairs):.1f}%); community: {ag_c}/{len(pairs)} ({100*ag_c/len(pairs):.1f}%)")
        print(f"  together at {A} but split at {B}: {len(ssa - ssb)}; together at {B} but split at {A}: {len(ssb - ssa)}")
    for i in sorted(a_all | b_all):
        rows.append({"replicon_group": g, "isolate": i, f"kept_{A}": i in a_kept, f"kept_{B}": i in b_kept, f"type_{A}": a_typ.get(i, ""), f"type_{B}": b_typ.get(i, "")})
out = f"pling_depth_comparison_{A}_vs_{B}.tsv"
if rows:
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rows)
    print(f"written: {out}")
