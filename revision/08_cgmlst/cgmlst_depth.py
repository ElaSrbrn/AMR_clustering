#!/usr/bin/env python3
"""Depth-dependent cgMLST: compare chewBBACA profiles at each depth with the reference depth.

Reference = profile of the same isolate at the highest available depth >= --min-ref-depth.
Called locus = EXC or INF allele (INF-<n> is treated as allele <n>); every other class
(LNF, PLOT3, PLOT5, LOTSC, NIPH, NIPHEM, ASM, ALM, PAMA, PLNF) counts as missing.

Per isolate and depth (cgmlst_vs_reference.tsv):
  called / missing loci by class, loci called in both, allele differences among those,
  loci lost (called in the reference, missing at this depth).
Pairwise (cgmlst_pairwise.tsv): allele distance between all isolate pairs of a scheme, counting
  loci called in both genomes, at each depth and at the reference; change in distance; with
  --thresholds, whether each pair is linked (distance <= threshold) at depth and reference.
Summary per scheme and depth, and pooled over schemes: cgmlst_summary.md

Usage:
  python3 cgmlst_depth.py --work cgmlst_work [--min-ref-depth 60] [--thresholds kpsc=15,ecoli=10]
"""
import argparse, csv, re, sys, statistics, collections, itertools
from pathlib import Path

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--work", default="cgmlst_work")
ap.add_argument("--min-ref-depth", type=int, default=60)
ap.add_argument("--thresholds", default="", help="scheme=alleles,... for cluster concordance (optional)")
ap.add_argument("--exclude", default="", help="comma-separated isolate numbers to leave out (sensitivity analysis)")
ap.add_argument("--out", default="cgmlst")
args = ap.parse_args()
work = Path(args.work)
thr = {k: int(v) for k, v in (x.split("=") for x in args.thresholds.split(",") if x)}

man = list(csv.DictReader(open(work / "manifest.tsv"), delimiter="\t"))
excl = {x.strip() for x in args.exclude.split(",") if x.strip()}
meta = {r["genome"]: r for r in man if r["id"] not in excl}
MISSING = ["LNF", "PLOT3", "PLOT5", "LOTSC", "NIPH", "NIPHEM", "ASM", "ALM", "PAMA", "PLNF"]

def parse(v):
    v = v.strip()
    m = re.fullmatch(r"(?:INF-)?\*?(\d+)", v)
    return ("call", m.group(1)) if m else ("miss", v.split("-")[0])

def q(xs):
    if not xs: return "NA"
    xs = sorted(xs)
    med = statistics.median(xs)
    lo, hi = (xs[int(0.25 * (len(xs) - 1))], xs[int(0.75 * (len(xs) - 1))])
    fmt = lambda x: f"{x:.1f}" if isinstance(x, float) and x != int(x) else f"{int(x)}"
    return f"{fmt(med)} ({fmt(lo)}\u2013{fmt(hi)}; max {fmt(xs[-1])})"

per_rows, pair_rows, summ = [], [], []
pooled = collections.defaultdict(lambda: collections.defaultdict(list))
for sdir in sorted(p for p in work.iterdir() if (p / "allelecall" / "results_alleles.tsv").exists()):
    scheme = sdir.name
    with open(sdir / "allelecall" / "results_alleles.tsv", newline="") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        loci = next(rdr)[1:]
        prof = {row[0]: [parse(v) for v in row[1:]] for row in rdr}
    nloc = len(loci)
    by_iso = collections.defaultdict(dict)                 # id -> depth -> genome
    for g in prof:
        m = meta.get(g)
        if m is None:
            if g.split("_")[0] not in excl:
                print(f"WARNING: {g} not in manifest", file=sys.stderr)
            continue
        by_iso[m["id"]][int(m["depth"])] = g
    ref = {}
    for i, ds in by_iso.items():
        ok = [d for d in ds if d >= args.min_ref_depth]
        if ok:
            ref[i] = ds[max(ok)]
    # per isolate x depth
    for i, ds in sorted(by_iso.items(), key=lambda x: int(x[0])):
        rg = ref.get(i)
        for d, g in sorted(ds.items(), reverse=True):
            p = prof[g]
            called = sum(t == "call" for t, _ in p)
            miss = collections.Counter(v for t, v in p if t == "miss")
            row = dict(scheme=scheme, id=i, s1_group=meta[g]["s1_group"], species=meta[g]["species"],
                       depth=d, genome=g, n_loci=nloc, called=called, pct_called=round(100 * called / nloc, 2),
                       **{f"n_{c}": miss.get(c, 0) for c in MISSING},
                       reference=rg or "none", shared="", allele_diff="", lost_vs_ref="", is_reference="no")
            if rg:
                pr = prof[rg]
                shared = [(a, b) for (ta, a), (tb, b) in zip(p, pr) if ta == tb == "call"]
                row.update(shared=len(shared), allele_diff=sum(a != b for a, b in shared),
                           lost_vs_ref=sum(ta != "call" and tb == "call" for (ta, _), (tb, _) in zip(p, pr)),
                           is_reference="yes" if g == rg else "no")
            per_rows.append(row)
    # pairwise distances, isolates with a reference
    ids = sorted(ref, key=int)
    def dist(g1, g2):
        a, b = prof[g1], prof[g2]
        sh = [(x, y) for (tx, x), (ty, y) in zip(a, b) if tx == ty == "call"]
        return sum(x != y for x, y in sh), len(sh)
    depths = sorted({d for i in ids for d in by_iso[i]}, reverse=True)
    for d in depths:
        for i1, i2 in itertools.combinations([i for i in ids if d in by_iso[i]], 2):
            dr, sr = dist(ref[i1], ref[i2])
            dd, sd = dist(by_iso[i1][d], by_iso[i2][d])
            row = dict(scheme=scheme, depth=d, id1=i1, id2=i2, dist_ref=dr, shared_ref=sr, dist_depth=dd,
                       shared_depth=sd, delta=dd - dr)
            if scheme in thr:
                row.update(threshold=thr[scheme], linked_ref=dr <= thr[scheme], linked_depth=dd <= thr[scheme])
            pair_rows.append(row)
    # summary per depth
    for d in depths:
        rs = [r for r in per_rows if r["scheme"] == scheme and r["depth"] == d]
        cmp_ = [r for r in rs if r["reference"] != "none" and r["is_reference"] == "no"]
        ps = [r for r in pair_rows if r["scheme"] == scheme and r["depth"] == d]
        line = dict(scheme=scheme, depth=d, genomes=len(rs), pct_called=q([r["pct_called"] for r in rs]),
                    compared=len(cmp_), diff=q([r["allele_diff"] for r in cmp_]),
                    zero=f"{sum(r['allele_diff'] == 0 for r in cmp_)}/{len(cmp_)}" if cmp_ else "NA",
                    lost=q([r["lost_vs_ref"] for r in cmp_]),
                    pairs=len(ps), pair_delta=q([abs(r["delta"]) for r in ps]) if ps else "NA")
        if scheme in thr and ps:
            lr = [r for r in ps if r["linked_ref"]]; ld = [r for r in ps if r["linked_depth"]]
            line["clusters"] = (f"linked at ref {len(lr)}, at depth {len(ld)}, both {sum(r['linked_depth'] for r in lr)}"
                                f" (threshold {thr[scheme]})")
        summ.append(line)
        for r in rs:
            pooled[d]["pct"].append(r["pct_called"])
        for r in cmp_:
            pooled[d]["diff"].append(r["allele_diff"]); pooled[d]["zero"].append(r["allele_diff"] == 0)

if not per_rows:
    sys.exit("ERROR: no results_alleles.tsv found under " + str(work))
for name, rows in ((f"{args.out}_vs_reference.tsv", per_rows), (f"{args.out}_pairwise.tsv", pair_rows)):
    if rows:
        keys = list(dict.fromkeys(k for r in rows for k in r))
        with open(name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t"); w.writeheader(); w.writerows(rows)

L = [f"## cgMLST by depth (reference: highest depth >= {args.min_ref_depth}x of the same isolate)", "",
     *( [f"Excluded isolates: {', '.join(sorted(excl, key=int))}", ""] if excl else [] ),
     "Called = EXC or INF. Values: median (IQR; max). Allele differences counted over loci called at both depths.", "",
     "| Scheme | Depth | Genomes | % loci called | Compared with reference | Allele differences | Isolates with 0 differences | Loci lost vs reference | Pairs | abs. change in pairwise distance |",
     "|---|---|---|---|---|---|---|---|---|---|"]
for s in summ:
    L.append(f"| {s['scheme']} | {s['depth']}x | {s['genomes']} | {s['pct_called']} | {s['compared']} | {s['diff']} | {s['zero']} | {s['lost']} | {s['pairs']} | {s['pair_delta']} |")
L += ["", "### All schemes pooled", "", "| Depth | Genomes | % loci called | Allele differences | Isolates with 0 differences |", "|---|---|---|---|---|"]
for d in sorted(pooled, reverse=True):
    z = pooled[d]["zero"]
    L.append(f"| {d}x | {len(pooled[d]['pct'])} | {q(pooled[d]['pct'])} | {q(pooled[d]['diff'])} | "
             f"{(str(sum(z)) + '/' + str(len(z))) if z else 'NA'} |")
cl = [s for s in summ if "clusters" in s]
if cl:
    L += ["", "### Cluster concordance", ""] + [f"- {s['scheme']} {s['depth']}x: {s['clusters']}" for s in cl]
noref = sorted({(r["scheme"], r["id"]) for r in per_rows if r["reference"] == "none"})
L += ["", f"Isolates without a reference >= {args.min_ref_depth}x (profiled, not compared): "
      + (", ".join(f"{i} ({s})" for s, i in noref) or "none")]
Path(f"{args.out}_summary.md").write_text("\n".join(L) + "\n")
print("\n".join(L))
print(f"\nwritten: {args.out}_summary.md, {args.out}_vs_reference.tsv, {args.out}_pairwise.tsv", file=sys.stderr)
