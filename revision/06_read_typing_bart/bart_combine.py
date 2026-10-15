#!/usr/bin/env python3
"""Collect bart MLST and AMR results from <reads_dir>/bart/<sample>/ into tables.

Usage:
  python3 bart_combine.py <out_prefix> <reads_dir> [<reads_dir> ...]
         [--catalog ReferenceGeneCatalog.txt] [--min-rel-depth 0.2]
         [--gold-dirs annotation_100x annotation_80x annotation_60x ...]
Writes <out_prefix>_mlst.tsv, <out_prefix>_amr.tsv and, with --gold-dirs, <out_prefix>_amr_vs_gold.tsv.

MLST: bart prints one line per sample, or several lines if more than one ST is equally
close (ambiguous). st_exact is 'yes' only if exactly one ST was reported and no allele
carries '?' (identity or coverage < 100), '-' (missing) or a differing top hit ')['.

AMR: every bart hit is kept; filter columns:
  rank_in_family  1 = highest-depth hit within its gene family in that sample. KMA assigns
                  single reads to neighbouring alleles (e.g. KPC-116 at ~1x next to KPC-2 at 8x);
                  these appear as rank >= 2.
  rel_depth       hit depth / median depth of the rank-1 hits of that sample.
  in_core_amr     with --catalog: type == AMR and scope == core, i.e. what AMRFinderPlus
                  reports without --plus (metal/biocide/virulence genes excluded).
  keep            rank_in_family == 1 and rel_depth >= --min-rel-depth (and in_core_amr if
                  a catalog is given). Does not use the gold standard.

Gold-standard comparison (Cottingham et al., Microbial Genomics 2026, doi:10.1099/mgen.0.001657):
for every gene family in the isolate's AMRFinderPlus assembly result, the top n bart hits by depth
are selected, where n is the number of DISTINCT alleles of that family in the assembly
(column keep_gold_topn). This measures allele concordance given known gene content.
--gold-dirs: annotation folders in priority order (highest depth first); for each isolate the
first folder containing it is used. Isolates are matched by sample name without the trailing
_<depth>x. Gold hits are restricted to Type AMR, excluding point mutations.
category: concordant | allele_mismatch | gold_not_exact (no exact allele call in the assembly)
| not_detected.

Carbapenemase split: the catalogue groups carbapenemase and non-carbapenemase alleles in some gene
families (e.g. blaOXA: OXA-48-like and OXA-1; blaGES: GES-5 and GES-1). In these mixed families,
carbapenemase alleles (subclass CARBAPENEM) form their own unit '<family>_carbapenemase' for ranking,
filtering and the gold-standard comparison, so that a co-carried non-carbapenemase allele cannot
mask a carbapenemase under the rank-1 rule. Requires --catalog.

Catalog annotation is joined on the reference nucleotide accession that bart prints at the
end of the description (e.g. 'NG 049253.1:6-887').
"""
import sys, re, csv, argparse, statistics, collections
from pathlib import Path

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("out_prefix")
ap.add_argument("dirs", nargs="+")
ap.add_argument("--catalog", help="ReferenceGeneCatalog.txt of the AMRFinderPlus database used for the index")
ap.add_argument("--min-rel-depth", type=float, default=0.2)
ap.add_argument("--gold-dirs", nargs="+", default=[],
                help="AMRFinderPlus annotation folders (annotation_<cov>), highest depth first")
args = ap.parse_args()

def allele_flagged(token):
    m = re.match(r"^[^(]+\((\d+)([?~-]?)", token)
    return m is None or bool(m.group(2)) or ")[" in token

ACC_RE = re.compile(r"\b([A-Z]{1,6})([ _]?)(\d+\.\d+):\d+-\d+\s*$")
def accession_from_description(desc):
    m = ACC_RE.search(desc)
    if not m:
        return ""
    prefix, sep, num = m.groups()
    return f"{prefix}_{num}" if sep else f"{prefix}{num}"

# ---- catalogue: nucleotide accession -> annotation ----
catalog, allele2fam, families = {}, {}, set()
fam_carb, fam_other = set(), set()
if args.catalog:
    with open(args.catalog, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        acc_cols = [c for c in rdr.fieldnames if "nucleotide_accession" in c.lower()]
        if not acc_cols:
            sys.exit(f"ERROR: no *nucleotide_accession column in {args.catalog}")
        for row in rdr:
            ann = {k: row.get(k, "") for k in ("allele", "gene_family", "type", "subtype",
                                               "class", "subclass", "scope", "product_name")}
            for c in acc_cols:
                a = (row.get(c) or "").strip()
                if a:
                    catalog[a] = ann
            if ann["allele"]:
                allele2fam[ann["allele"]] = ann["gene_family"] or ann["allele"]
            if ann["gene_family"]:
                families.add(ann["gene_family"])
                (fam_carb if "CARBAPENEM" in ann["subclass"].upper() else fam_other).add(ann["gene_family"])
    print(f"catalog: {len(catalog)} accessions", file=sys.stderr)

MIXED = fam_carb & fam_other
if MIXED:
    print(f"families split into carbapenemase / other alleles: {len(MIXED)}", file=sys.stderr)
def fam_key(family, subclass):
    """ranking unit: gene family, with carbapenemase alleles of mixed families as their own unit"""
    if family in MIXED and "CARBAPENEM" in (subclass or "").upper():
        return family + "_carbapenemase"
    return family

DEPTH_SUFFIX = re.compile(r"_\d+x$")
def isolate_key(sample):
    return DEPTH_SUFFIX.sub("", sample)

def symbol_family(sym, subclass=""):
    return fam_key(allele2fam.get(sym, sym), subclass)

# ---- gold standard: AMRFinderPlus on assemblies ----
gold = {}
for gdir in args.gold_dirs:
    gdir = Path(gdir).resolve()
    if not gdir.is_dir():
        sys.exit(f"ERROR: gold dir not found: {gdir}")
    for f in sorted(gdir.glob("*/amrfinder/*_amrfinder.tsv")):
        key = isolate_key(f.parent.parent.name)
        if key in gold:
            continue   # a higher-depth folder already provided this isolate
        fams = {}
        with f.open(newline="") as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            fn = rdr.fieldnames or []
            col = lambda *names: next((c for c in fn for n in names if c.lower() == n.lower()), None)
            sym_c, type_c, sub_c = col("Element symbol", "Gene symbol"), col("Type", "Element type"), col("Subtype", "Element subtype")
            meth_c, id_c, cov_c = col("Method"), col("% Identity to reference"), col("% Coverage of reference")
            subcl_c = col("Subclass")
            for row in rdr:
                if type_c and row.get(type_c) != "AMR":
                    continue
                if sub_c and row.get(sub_c) == "POINT":
                    continue
                sym = (row.get(sym_c) or "").strip()
                if sym:
                    fams.setdefault(symbol_family(sym, row.get(subcl_c, "") if subcl_c else ""), []).append(dict(
                        sym=sym, method=row.get(meth_c, "") if meth_c else "",
                        ident=row.get(id_c, "") if id_c else "", cov=row.get(cov_c, "") if cov_c else ""))
        gold[key] = dict(source=gdir.name, fams=fams)
if args.gold_dirs:
    print(f"gold standard: {len(gold)} isolates", file=sys.stderr)

gold_rows, no_gold = [], set()
mlst_rows, amr_rows, missing, unmatched = [], [], [], set()
for d in args.dirs:
    d = Path(d).resolve()
    group = f"{d.parent.name}/{d.name}"
    for sdir in sorted((d / "bart").glob("*/")):
        sample = sdir.name
        mf, af = sdir / f"{sample}_mlst.tsv", sdir / f"{sample}_amr.tsv"

        # ---- MLST ----
        lines = [l.split("\t") for l in mf.read_text().splitlines() if l.strip()] if mf.exists() else []
        if not lines:
            missing.append(f"{group}\t{sample}\tmlst")
            mlst_rows.append(dict(group=group, sample=sample, scheme="", st="", n_st_reported=0,
                                  n_imperfect_alleles="", st_exact="no", alleles=""))
        else:
            sts = sorted({l[2] for l in lines if len(l) > 2})
            first = lines[0]
            alleles = [a for a in first[3:] if not re.match(r"(clonal_complex|species|CC|Lineage)\(", a)]
            imperfect = sum(1 for a in alleles if allele_flagged(a))
            exact = "yes" if len(sts) == 1 and sts[0] != "?" and imperfect == 0 else "no"
            mlst_rows.append(dict(group=group, sample=sample, scheme=first[1] if len(first) > 1 else "",
                                  st=";".join(sts), n_st_reported=len(sts), n_imperfect_alleles=imperfect,
                                  st_exact=exact, alleles=" ".join(alleles)))

        # ---- AMR ----
        if not af.exists():
            missing.append(f"{group}\t{sample}\tamr")
            continue
        hits = []
        for r in list(csv.reader(af.open(), delimiter="\t"))[1:]:
            if len(r) < 7:
                continue
            acc = accession_from_description(r[2])
            ann = catalog.get(acc, {})
            if args.catalog and not ann:
                unmatched.add(acc or r[2][:60])
            hits.append(dict(group=group, sample=sample, gene=r[1], description=r[2], accession=acc,
                             allele=ann.get("allele", ""),
                             gene_family=fam_key(ann.get("gene_family", "") or r[1], ann.get("subclass", "")),
                             type=ann.get("type", ""), scope=ann.get("scope", ""), amr_class=ann.get("class", ""),
                             subclass=ann.get("subclass", ""), length=r[3], identity=float(r[4]),
                             coverage=float(r[5]), depth=float(r[6])))
        if not hits:
            amr_rows.append(dict(group=group, sample=sample, gene="", description="NO_HITS"))
            continue

        fam = {}
        for h in sorted(hits, key=lambda x: (-x["depth"], -x["identity"])):
            fam.setdefault(h["gene_family"], []).append(h)
        for members in fam.values():
            for i, h in enumerate(members, 1):
                h["rank_in_family"] = i
        ref = statistics.median([m[0]["depth"] for m in fam.values()])
        for h in hits:
            h["rel_depth"] = round(h["depth"] / ref, 3) if ref > 0 else ""
            core = (h["type"] == "AMR" and h["scope"] == "core") if args.catalog else True
            h["in_core_amr"] = ("yes" if core else "no") if args.catalog else ""
            h["keep"] = "yes" if (h["rank_in_family"] == 1 and ref > 0 and
                                  h["depth"] / ref >= args.min_rel_depth and core) else "no"

        # ---- comparison with assembly gold standard ----
        if args.gold_dirs:
            g = gold.get(isolate_key(sample))
            for h in hits:
                h["keep_gold_topn"] = "no" if g else "NA"
            if not g:
                no_gold.add(isolate_key(sample))
            else:
                for gfam, ghits in sorted(g["fams"].items()):
                    galleles = sorted({x["sym"] for x in ghits})      # distinct alleles; copies collapse
                    n = len(galleles)
                    top = fam.get(gfam, [])[:n]
                    for h in top:
                        h["keep_gold_topn"] = "yes"
                    bart_alleles = sorted({h["allele"] or h["gene"] for h in top})
                    methods = sorted({x["method"] for x in ghits})
                    gold_exact = all(methods) and all(m.startswith(("EXACT", "ALLELE")) for m in methods)
                    if not top:
                        category = "not_detected"
                    elif galleles == bart_alleles:
                        category = "concordant"
                    elif not gold_exact:
                        category = "gold_not_exact"
                    else:
                        category = "allele_mismatch"
                    gold_rows.append(dict(
                        group=group, sample=sample, isolate=isolate_key(sample), gold_source=g["source"],
                        gene_family=gfam, n_gold=n, n_gold_copies=len(ghits),
                        gold_alleles=";".join(galleles), gold_methods=";".join(methods),
                        gold_identity=";".join(x["ident"] for x in ghits),
                        gold_coverage=";".join(x["cov"] for x in ghits),
                        bart_alleles=";".join(bart_alleles),
                        bart_depths=";".join(f'{h["depth"]:.2f}' for h in top),
                        bart_identities=";".join(f'{h["identity"]:.2f}' for h in top),
                        family_detected="yes" if top else "no",
                        allele_concordant="yes" if galleles == bart_alleles else "no",
                        category=category))
        amr_rows.extend(hits)

# ---- write tables ----
def write(path, rows, fields=None):
    if not rows:
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields or list(rows[0].keys()), delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"{path}: {len(rows)} rows", file=sys.stderr)

write(f"{args.out_prefix}_mlst.tsv", mlst_rows)
write(f"{args.out_prefix}_amr.tsv", amr_rows,
      ["group", "sample", "gene", "gene_family", "allele", "type", "scope", "amr_class", "subclass",
       "identity", "coverage", "depth", "rel_depth", "rank_in_family", "in_core_amr", "keep",
       "keep_gold_topn", "accession", "length", "description"])
if args.gold_dirs:
    write(f"{args.out_prefix}_amr_vs_gold.tsv", gold_rows)
    cats = collections.Counter(r["category"] for r in gold_rows)
    print("categories: " + ", ".join(f"{k} {v}" for k, v in sorted(cats.items())), file=sys.stderr)
    if no_gold:
        print(f"WARNING: {len(no_gold)} isolates without gold standard: {sorted(no_gold)}", file=sys.stderr)
if unmatched:
    print(f"WARNING: {len(unmatched)} accessions not in catalog, e.g. {sorted(unmatched)[:5]}", file=sys.stderr)
if missing:
    with open(f"{args.out_prefix}_missing.tsv", "w") as fh:
        fh.write("group\tsample\toutput\n" + "\n".join(missing) + "\n")
    print(f"{len(missing)} missing outputs, listed in {args.out_prefix}_missing.tsv", file=sys.stderr)
