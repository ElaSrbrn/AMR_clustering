#!/usr/bin/env python3
"""Summarise read-based (bart) typing and compare it with assembly-based typing at the same depth.

Both methods are scored against the same gold standard: AMRFinderPlus (and mlst) on the assembly
of the same isolate at the highest available depth (>= 40x; list --gold-dirs highest depth first,
the first folder containing an isolate is used). Only isolates present in both the read-based
and the assembly-based results at a given depth are compared (paired design).

AMR (per isolate x gold gene family; core AMR genes, point mutations excluded):
  family detected     the method reports at least one allele of the gold family
  allele concordant   the set of distinct alleles equals the gold set
  rows where the gold standard has no exact allele call (AMRFinderPlus method not EXACT/ALLELE)
  are excluded from allele concordance and counted separately.
  Read-based calls use the top-n rule (bart_combine.py, column category).
Carbapenemase genes (gold allele with subclass CARBAPENEM in the catalogue) are summarised separately; the isolate-level outcome "all carbapenemase
alleles correct" is compared between methods with an exact McNemar test.

False positives (carbapenemases, with --reads-hits): calls are made WITHOUT the gold standard.
  reads, unfiltered     every carbapenemase hit counts
  reads, keep >= t      highest-depth hit per gene family, core AMR, depth >= t x median depth
                        of the isolate's gene families (t from --thresholds)
  assembly              every carbapenemase reported by AMRFinderPlus at that depth
  Family level: FP = carbapenemase family called but absent from the gold standard.
  Secondary alleles: rank >= 2 carbapenemase hits within a gold family whose allele is not in the
  gold standard (e.g. KPC-116 next to KPC-2), read rules only. The keep rule reports rank 1 only,
  so these show what a depth threshold without the rank restriction would add.
  Isolate level, gold-independent: "all carbapenemase alleles correct" for reads under
  keep >= --keep-threshold vs assembly, exact McNemar test.
  Isolate level: sensitivity and specificity for "carbapenemase present".
  Gene families of the catalogue that mix carbapenemase and other alleles (e.g. blaOXA, blaGES)
  are split as in bart_combine.py ('<family>_carbapenemase').
  Details of every FP are written to <out>_carbapenemase_fp.tsv; carbapenemase calls of either
  method that differ from the gold standard to <out>_carbapenemase_discordant.tsv.

MLST: ST from bart (exact calls only) and from mlst on the assembly at the same depth, each
compared with the mlst ST of the highest-depth assembly with a numeric ST and depth >= --min-gold-depth.

Usage:
  python3 compare_reads_vs_assemblies.py \
      --reads-vs-gold bart_results_amr_vs_gold.tsv --bart-mlst bart_results_mlst.tsv \
      --catalog ReferenceGeneCatalog.txt \
      --gold-dirs carblot_pos/annotation_100x carblot_pos/annotation_80x carblot_pos/annotation_60x \
                  carblot_pos/annotation_40x \
                  carblot_neg/annotation_100x carblot_neg/annotation_60x carblot_neg/annotation_40x \
      [--reads-hits bart_results_amr.tsv] [--thresholds 0.1,0.2,0.3,0.5] \
      [--assembly-mlst all_coverages_mlst.tsv] [--min-gold-depth 40] [--out summary]
With --isolate-groups (from classify_isolates.py), results are reported per manuscript group
(e.g. S1a/S1b) instead of per sequencing folder, and isolates not listed are excluded.
Assembly annotations at the read depths are expected at <set>/annotation_<depth> (e.g.
carblot_pos/annotation_10x for the read group carblot_pos/10x).
"""
import argparse, csv, re, sys, math, collections
from pathlib import Path

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--reads-vs-gold", required=True)
ap.add_argument("--bart-mlst")
ap.add_argument("--catalog", required=True)
ap.add_argument("--gold-dirs", nargs="+", required=True)
ap.add_argument("--reads-hits", help="bart_results_amr.tsv (all read hits) for the false-positive analysis")
ap.add_argument("--thresholds", default="0.1,0.2,0.3,0.5",
                help="relative-depth thresholds for the gold-independent read filter")
ap.add_argument("--keep-threshold", type=float, default=0.2,
                help="threshold of the gold-independent read rule used for the isolate-level McNemar test")
ap.add_argument("--assembly-mlst", help="mlst table with columns coverage, isolate, mlst_st (07_mlst.sh output)")
ap.add_argument("--min-gold-depth", type=int, default=40,
                help="lowest assembly depth accepted as MLST reference (x); use the same floor as --gold-dirs")
ap.add_argument("--isolate-groups",
                help="tsv isolate<TAB>group (classify_isolates.py); results are grouped by these groups "
                     "instead of the sequencing folders, isolates not listed are excluded")
ap.add_argument("--out", default="summary")
args = ap.parse_args()

CARB = re.compile(r"^bla(KPC|NDM|VIM|IMP|OXA-48|OXA-23|OXA-24|OXA-58|GES|IMI|SME|FRI|GIM|SIM|DIM|SPM|AIM)")
DEPTH = re.compile(r"_(\d+x)$")
key = lambda s: DEPTH.sub("", s)

allele2fam, allele2sub = {}, {}
fam_carb, fam_other = set(), set()
with open(args.catalog, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r.get("allele"):
            allele2fam[r["allele"]] = r.get("gene_family") or r["allele"]
            allele2sub[r["allele"]] = (r.get("subclass") or "").upper()
        if r.get("gene_family"):
            (fam_carb if "CARBAPENEM" in (r.get("subclass") or "").upper() else fam_other).add(r["gene_family"])
MIXED = fam_carb & fam_other
SPLIT = "_carbapenemase"

def fam_key(family, subclass):
    """same rule as bart_combine.py: carbapenemase alleles of mixed families form their own unit"""
    if family in MIXED and "CARBAPENEM" in (subclass or "").upper():
        return family + SPLIT
    return family

def is_carb(row):
    """Carbapenemase if any gold allele has subclass CARBAPENEM in the AMRFinderPlus catalogue
    (so e.g. ESBL variants of GES are not counted); family-name pattern only as fallback."""
    if row["gene_family"].endswith(SPLIT):
        return True
    alleles = [a for a in row.get("gold_alleles", "").split(";") if a]
    known = [a for a in alleles if a in allele2sub]
    if known:
        return any("CARBAPENEM" in allele2sub[a] for a in known)
    return row["gene_family"] not in MIXED and bool(CARB.match(row["gene_family"]))
fam_of = lambda s: allele2fam.get(s, s)

def is_carb_symbol(sym, subclass=""):
    if subclass:
        return "CARBAPENEM" in subclass.upper()
    if sym in allele2sub:
        return "CARBAPENEM" in allele2sub[sym]
    return bool(CARB.match(fam_of(sym)))

def read_amrfinder(f):
    fams = collections.defaultdict(lambda: {"alleles": set(), "methods": set(), "carb": False})
    with open(f, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        fn = rdr.fieldnames or []
        col = lambda *n: next((c for c in fn for x in n if c.lower() == x.lower()), None)
        sc, tc, uc, mc = col("Element symbol", "Gene symbol"), col("Type", "Element type"), col("Subtype", "Element subtype"), col("Method")
        ssc = col("Subclass")
        for r in rdr:
            if (tc and r.get(tc) != "AMR") or (uc and r.get(uc) == "POINT"):
                continue
            s = (r.get(sc) or "").strip()
            if s:
                sub = (r.get(ssc) or "") if ssc else allele2sub.get(s, "")
                f = fam_key(fam_of(s), sub)
                fams[f]["alleles"].add(s)
                fams[f]["methods"].add(r.get(mc, "") if mc else "")
                if ("CARBAPENEM" in sub.upper()) if ssc else is_carb_symbol(s):
                    fams[f]["carb"] = True
    return fams

def load_dir(d):
    out = {}
    for f in sorted(Path(d).glob("*/amrfinder/*_amrfinder.tsv")):
        out[key(f.parent.parent.name)] = read_amrfinder(f)
    return out

# ---- optional regrouping by manuscript classification instead of sequencing folder ----
ISO_GROUP = {}
if args.isolate_groups:
    with open(args.isolate_groups, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ISO_GROUP[r["isolate"]] = r["group"]
FOLDER_SETS = set()

def regroup(r, sample_field):
    """returns False if the row is to be dropped; otherwise sets r['group'] to <class>/<depth>"""
    sset, depth = r["group"].split("/")
    FOLDER_SETS.add(sset)
    if not ISO_GROUP:
        return True
    g = ISO_GROUP.get(key(r[sample_field]))
    if g is None:
        return False
    r["group"] = f"{g}/{depth}"
    return True

def annotation_dirs(group):
    sset, depth = group.split("/")
    sets = sorted(FOLDER_SETS) if ISO_GROUP else [sset]
    return [Path(s) / f"annotation_{depth}" for s in sets if (Path(s) / f"annotation_{depth}").is_dir()]

def load_dirs(dirs):
    out = {}
    for d in dirs:
        for k, v in load_dir(d).items():
            out.setdefault(k, v)
    return out

gold = {}
for d in args.gold_dirs:
    for k, v in load_dir(d).items():
        gold.setdefault(k, v)
exact = lambda methods: all(methods) and all(m.startswith(("EXACT", "ALLELE")) for m in methods)

# ---- read-based rows (already scored by bart_combine.py) ----
reads = collections.defaultdict(list)          # group -> rows
with open(args.reads_vs_gold, newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
if "category" not in rows[0]:
    sys.exit("ERROR: reads-vs-gold table has no 'category' column; rerun the current bart_combine.py")
for r in rows:
    r["gx"] = exact([m for m in r.get("gold_methods", "").split(";")] or [""])
    if regroup(r, "sample"):
        reads[r["group"]].append(r)
if ISO_GROUP:
    print(f"grouped by --isolate-groups: {len({key(r['sample']) for v in reads.values() for r in v})} isolates "
          f"in groups {sorted(reads)}", file=sys.stderr)

# ---- assembly-based rows at the same depth, scored the same way ----
def score_assembly(group):
    dirs = annotation_dirs(group)
    if not dirs:
        return None
    asm = load_dirs(dirs)
    out = []
    for iso, fams in asm.items():
        g = gold.get(iso)
        if not g:
            continue
        for gfam, gv in g.items():
            a = fams.get(gfam, {"alleles": set()})["alleles"]
            if not a:
                cat = "not_detected"
            elif a == gv["alleles"]:
                cat = "concordant"
            elif not exact(gv["methods"]):
                cat = "gold_not_exact"
            else:
                cat = "allele_mismatch"
            out.append(dict(isolate=iso, gene_family=gfam, category=cat, gx=exact(gv["methods"]),
                            called=";".join(sorted(a)), methods=";".join(sorted(fams.get(gfam, {}).get("methods", set()))),
                            gold_alleles=";".join(sorted(gv["alleles"]))))
    return out

def metrics(rs):
    """family detection over all gold families; allele concordance over detected families whose
    gold call is exact (same rows for reads and assembly)"""
    n = len(rs)
    det = sum(r["category"] != "not_detected" for r in rs)
    evaluable = [r for r in rs if r["category"] != "not_detected" and r["gx"]]
    conc = sum(r["category"] == "concordant" for r in evaluable)
    pct = lambda a, b: f"{a}/{b} ({100*a/b:.1f}%)" if b else "NA"
    return pct(det, n), pct(conc, len(evaluable)), sum(not r["gx"] for r in rs)

def isolate_carb_ok(rs):
    """isolate -> True if every gold carbapenemase family is detected with the correct allele set."""
    ok = {}
    for r in rs:
        if is_carb(r):
            ok.setdefault(r["isolate"], True)
            if r["category"] == "not_detected" or (r["gx"] and r["category"] != "concordant"):
                ok[r["isolate"]] = False
    return ok

def mcnemar_exact(b, c):
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    return min(1.0, 2 * sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n)

lines = ["## AMR: read-based (bart) vs assembly-based (AMRFinderPlus) at the same depth",
         "Paired isolates only. Reads: top-n rule (n from the gold standard). Allele concordance: detected "
         "families whose gold call is exact (EXACT/ALLELE); Gold not exact = gold families without such a call.",
         "",
         "| Group | Method | Isolates | Carbapenemase family detected | Carbapenemase allele concordant | All core AMR family detected | All core AMR allele concordant | Gold not exact |",
         "|---|---|---|---|---|---|---|---|"]
tests, asm_ok, disc = [], {}, []
for group in sorted(reads):
    rr = [dict(r, isolate=r["isolate"]) for r in reads[group]]
    ar = score_assembly(group)
    if ar is None:
        lines.append(f"| {group} | (no assembly annotation folder found) | | | | | | |")
        continue
    paired = {r["isolate"] for r in rr} & {r["isolate"] for r in ar}
    rr = [r for r in rr if r["isolate"] in paired]
    ar = [r for r in ar if r["isolate"] in paired]
    for name, rs in (("reads (bart)", rr), ("assembly", ar)):
        cd, cc, _ = metrics([r for r in rs if is_carb(r)])
        ad, ac, gne = metrics(rs)
        lines.append(f"| {group} | {name} | {len(paired)} | {cd} | {cc} | {ad} | {ac} | {gne} |")
    for name, rs in (("reads (bart)", rr), ("assembly", ar)):
        for r in rs:
            if is_carb(r) and r["category"] != "concordant":
                disc.append(dict(group=group, isolate=r["isolate"], method=name, gene_family=r["gene_family"],
                                 category=r["category"], gold_alleles=r["gold_alleles"],
                                 called=r.get("called", r.get("bart_alleles", "")),
                                 detail=r.get("methods") or r.get("bart_depths", "")))
    ro, ao = isolate_carb_ok(rr), isolate_carb_ok(ar)
    asm_ok[group] = ao
    both = set(ro) & set(ao)
    b = sum(ro[i] and not ao[i] for i in both)      # reads correct, assembly not
    c = sum(ao[i] and not ro[i] for i in both)      # assembly correct, reads not
    tests.append(f"| {group} | {len(both)} | {sum(ro[i] for i in both)} | {sum(ao[i] for i in both)} | {b} | {c} | {mcnemar_exact(b, c):.3g} |")

if disc:
    with open(f"{args.out}_carbapenemase_discordant.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(disc[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(disc)
    lines.append(f"\nCarbapenemase calls that differ from the gold standard: {args.out}_carbapenemase_discordant.tsv "
                 "(detail = AMRFinderPlus method for assemblies, read depths for bart)")
lines += ["", "### Isolate level: all carbapenemase alleles correct",
          "| Group | Isolates with carbapenemase | Reads correct | Assembly correct | Reads only | Assembly only | McNemar exact p |",
          "|---|---|---|---|---|---|---|"] + tests

# ---- false positives: carbapenemases, calls made without the gold standard ----
gold_carb = {iso: {f for f, v in fams.items() if v["carb"]} for iso, fams in gold.items()}

if args.reads_hits:
    thresholds = [float(t) for t in args.thresholds.split(",") if t]
    hits = collections.defaultdict(lambda: collections.defaultdict(list))   # group -> isolate -> hits
    with open(args.reads_hits, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not r.get("gene_family") or not r.get("depth"):
                continue
            if not regroup(r, "sample"):
                continue
            if r["gene_family"].endswith(SPLIT) or is_carb_symbol(r.get("allele") or r["gene"], r.get("subclass", "")):
                hits[r["group"]][key(r["sample"])].append(r)

    def read_calls(isolate_hits, t):
        """families called and extra alleles (rank >= 2) under rule t (None = unfiltered)"""
        fams, extra = set(), []
        for h in isolate_hits:
            rel = float(h["rel_depth"]) if h.get("rel_depth") not in (None, "") else 0.0
            rank = int(h.get("rank_in_family") or 1)
            passes = t is None or (rel >= t and h.get("in_core_amr", "yes") == "yes")
            if passes and (t is None or rank == 1):
                fams.add(h["gene_family"])
            if passes and rank >= 2:
                extra.append(h)
        return fams, extra

    fp_rows = []
    lines += ["", "## False positives: carbapenemases (calls made without the gold standard)",
              "Paired isolates with a gold standard. Family FP = carbapenemase family called but absent "
              "from the gold standard. Secondary alleles = hits of rank >= 2 within a gold family, allele not in "
              "the gold standard, passing the depth threshold (the keep rule itself reports only rank 1, so these "
              "show what a depth filter without the rank restriction would add).",
              "", "| Group | Method | Isolates | Gold carbapenemase families | Detected | Missed | FP families | Secondary alleles | Sensitivity (isolate) | Specificity (isolate) |",
              "|---|---|---|---|---|---|---|---|---|---|"]
    for group in sorted(reads):
        asm = load_dirs(annotation_dirs(group))
        isolates = sorted(i for i in ({key(r["sample"]) for r in reads[group]} & set(asm)) if i in gold_carb)
        if not isolates:
            continue
        rules = [("reads, unfiltered", None)] + [(f"reads, keep >= {t}", t) for t in thresholds] + [("assembly", "asm")]
        for name, t in rules:
            tp = fn = fp = extra_n = 0
            sens_hit = sens_n = spec_hit = spec_n = 0
            for iso in isolates:
                g = gold_carb[iso]
                if t == "asm":
                    called = {f for f, v in asm[iso].items() if v["carb"]}
                    extra = []
                else:
                    called, extra = read_calls(hits[group].get(iso, []), t)
                    extra = [h for h in extra if h["gene_family"] in g and
                             (h.get("allele") or h["gene"]) not in gold[iso][h["gene_family"]]["alleles"]]
                tp += len(called & g); fn += len(g - called); fp += len(called - g); extra_n += len(extra)
                for f in sorted(called - g):
                    det = [h for h in hits[group].get(iso, []) if h["gene_family"] == f] if t != "asm" else []
                    fp_rows.append(dict(group=group, isolate=iso, method=name, type="family_FP", gene_family=f,
                                        alleles=";".join(sorted({h.get("allele") or h["gene"] for h in det})) if det
                                        else ";".join(sorted(asm[iso][f]["alleles"])) if t == "asm" else "",
                                        depth=";".join(h["depth"] for h in det), rel_depth=";".join(h.get("rel_depth", "") for h in det)))
                for h in extra:
                    fp_rows.append(dict(group=group, isolate=iso, method=name, type="extra_allele", gene_family=h["gene_family"],
                                        alleles=h.get("allele") or h["gene"], depth=h["depth"], rel_depth=h.get("rel_depth", "")))
                if g:
                    sens_n += 1; sens_hit += bool(called)
                else:
                    spec_n += 1; spec_hit += not called
            pct = lambda a, b: f"{a}/{b} ({100*a/b:.1f}%)" if b else "NA"
            lines.append(f"| {group} | {name} | {len(isolates)} | {tp+fn} | {tp} | {fn} | {fp} | "
                         f"{extra_n if t != 'asm' else 'NA'} | {pct(sens_hit, sens_n)} | {pct(spec_hit, spec_n)} |")
    # isolate level, reads called WITHOUT the gold standard (keep rule) vs assembly
    K = args.keep_threshold
    lines += ["", f"### Isolate level: all carbapenemase alleles correct, reads with keep >= {K} (no gold information)",
              "Correct = every gold carbapenemase family called (rank 1, depth >= threshold, core AMR) with an allele "
              "contained in the gold allele set (presence only for gold calls that are not exact).",
              "| Group | Isolates with carbapenemase | Reads correct | Assembly correct | Reads only | Assembly only | McNemar exact p |",
              "|---|---|---|---|---|---|---|"]
    for group in sorted(asm_ok):
        ao = asm_ok[group]
        ko = {}
        for iso in ao:
            ok = True
            for f in gold_carb.get(iso, ()):
                gv = gold[iso][f]
                call = [h for h in hits[group].get(iso, []) if h["gene_family"] == f
                        and int(h.get("rank_in_family") or 1) == 1 and h.get("in_core_amr", "yes") == "yes"
                        and h.get("rel_depth") not in (None, "") and float(h["rel_depth"]) >= K]
                if not call or (exact(gv["methods"]) and (call[0].get("allele") or call[0]["gene"]) not in gv["alleles"]):
                    ok = False
            ko[iso] = ok
        b = sum(ko[i] and not ao[i] for i in ao); c = sum(ao[i] and not ko[i] for i in ao)
        lines.append(f"| {group} | {len(ao)} | {sum(ko.values())} | {sum(ao.values())} | {b} | {c} | {mcnemar_exact(b, c):.3g} |")

    if fp_rows:
        with open(f"{args.out}_carbapenemase_fp.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(fp_rows[0].keys()), delimiter="\t")
            w.writeheader(); w.writerows(fp_rows)
        lines.append(f"\nDetails of all false positives: {args.out}_carbapenemase_fp.tsv")

# ---- MLST ----
if args.bart_mlst:
    with open(args.bart_mlst, newline="") as fh:
        bm = [r for r in csv.DictReader(fh, delimiter="\t") if regroup(r, "sample")]
    asm_st, conflicts = {}, []
    if args.assembly_mlst:
        # the mlst table can hold repeated rows per isolate and depth; take the most frequent
        # non-empty ST and report isolates whose repeated rows disagree
        seen = collections.defaultdict(list)
        with open(args.assembly_mlst, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                st = (r.get("mlst_st") or "").strip()
                if st:
                    seen[(key(r["isolate"]), r.get("coverage", ""))].append(st)
        for k, sts in seen.items():
            c = collections.Counter(sts)
            asm_st[k] = c.most_common(1)[0][0]
            if len(c) > 1:
                conflicts.append(f"{k[0]} {k[1]}: {dict(c)}")
        if conflicts:
            print(f"WARNING: {len(conflicts)} isolate/depth combinations with conflicting mlst rows, e.g. {conflicts[:3]}",
                  file=sys.stderr)
    gold_st = {}
    for (iso, cov), st in asm_st.items():
        d = int(cov.rstrip("x")) if cov.rstrip("x").isdigit() else 0
        if st.isdigit() and d >= args.min_gold_depth and d >= gold_st.get(iso, (0, ""))[0]:
            gold_st[iso] = (d, st)
    gd = collections.Counter(f"{d}x" for d, _ in gold_st.values())
    lines += ["", "## MLST",
              "Exact = single ST with all alleles 100% identity and coverage (bart) or a numeric ST (mlst).",
              f"Reference ST: mlst on the highest-depth assembly >= {args.min_gold_depth}x with a numeric ST "
              f"({', '.join(f'{k}: {v}' for k, v in sorted(gd.items(), key=lambda x: -int(x[0][:-1])))} isolates).",
              "", "| Group | Samples | No scheme | Reads exact ST | Reads ST = gold | Assembly ST = gold |",
              "|---|---|---|---|---|---|"]
    by = collections.defaultdict(list)
    for r in bm:
        by[r["group"]].append(r)
    for group, rs in sorted(by.items()):
        depth = group.split("/")[1]
        typed = [r for r in rs if r["scheme"]]
        ex = [r for r in typed if r["st_exact"] == "yes"]
        with_gold = [r for r in typed if key(r["sample"]) in gold_st]
        rg = sum(r["st_exact"] == "yes" and r["st"] == gold_st[key(r["sample"])][1] for r in with_gold)
        ag = sum(asm_st.get((key(r["sample"]), depth), "") == gold_st[key(r["sample"])][1] for r in with_gold)
        pct = lambda a, b: f"{a}/{b} ({100*a/b:.1f}%)" if b else "NA"
        lines.append(f"| {group} | {len(rs)} | {len(rs)-len(typed)} | {pct(len(ex), len(typed))} | "
                     f"{pct(rg, len(with_gold))} | {pct(ag, len(with_gold)) if asm_st else 'NA'} |")

Path(f"{args.out}.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
print(f"\nwritten: {args.out}.md", file=sys.stderr)
