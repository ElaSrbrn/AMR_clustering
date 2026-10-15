#!/usr/bin/env python3
"""Match the isolates of the bart analysis to the manuscript classification (Supplementary Table S1).

S1a: carbapenemase-producing isolates (RESIST-5 positive, n = 100)
S1b: exploratory cohort, carbapenem-non-susceptible and RESIST-5 negative (n = 10)
The analysis folders carblot_pos/ and carblot_neg/ are sequencing batches and are not used as groups.

Writes
  isolate_classification.tsv   every analysed isolate: folder, S1 group, reference, carbapenemase genes
                               in the reference vs. S1 ('WGS' column), duplicates
  isolate_groups.tsv           isolate -> S1a/S1b, one run per isolate (input for
                               compare_reads_vs_assemblies.py --isolate-groups)
Usage (in the folder with carblot_pos/ and carblot_neg/):
  python3 classify_isolates.py --catalog ReferenceGeneCatalog.txt [--vs-gold bart_results_amr_vs_gold.tsv]
"""
import argparse, csv, re, sys, collections
from pathlib import Path

S1A = {
    "3989": ["blaKPC-2", "blaOXA-48"], "4704": ["blaVIM-1 (not detected at 10x)"],
    "4838": ["blaNDM-1"], "4920": ["blaOXA-48", "blaVIM-1"], "5043": ["NDM-1, blaDHA"],
    "5244": ["blaOXA-48"], "5305": ["blaVIM-86 (subtype not detected at 10x)"], "5399": ["blaKPC-2"],
    "5590": ["blaVIM-1"], "5860": ["blaKPC-2"], "5871": ["blaKPC-2"], "5969": ["blaNDM-1"],
    "5970": ["blaOXA-48", "blaNDM-1"], "5976": ["blaKPC-2"], "6009": ["blaKPC-3"],
    "6014": ["blaKPC-2"], "6035": ["blaNDM-1"], "6056": ["blaNDM-1"], "6057": ["blaNDM-5"],
    "6067": ["blaKPC-2"], "6068": ["blaKPC-2"], "6233": ["blaOXA-48"], "6243": ["blaKPC-2"],
    "6244": ["blaOXA-181 (subtype not detected at 10x)"], "6247": ["blaOXA-181"],
    "6287": ["blaOXA-48"], "6312": ["blaOXA-244"], "6365": ["blaKPC-2"], "6394": ["blaNDM-1"],
    "6378": ["blaNDM-5"], "6403": ["blaNDM-5"], "6433": ["blaOXA-181 (subtype not detected at 10x)"],
    "6455": ["blaOXA-181"], "6454": ["blaKPC-2 ( subtype not detected at 10x)"],
    "6472": ["blaOXA-48"], "6466": ["blaNDM-1"], "6482": ["blaKPC-2"], "6500": ["blaKPC-3"],
    "6520": ["blaVIM-1"], "6567": ["blaKPC-2"], "6571": ["blaOXA-48"], "6601": ["blaKPC-2"],
    "6623": ["blaVIM-1"], "6669": ["blaKPC-2"], "6670": ["blaKPC-2"], "6701": ["blaVIM-1"],
    "6711": ["blaKPC-2"], "6712": ["blaKPC-2"], "6726": ["blaKPC-2"], "6738": ["blaOXA-48"],
    "6735": ["blaOXA-48"], "6737": ["blaOXA-232", "blaNDM-5"], "6748": ["blaNDM-5"],
    "6779": ["blaOXA-48 (not detected at 10x)"], "6774": ["blaOXA-48"], "6775": ["blaOXA-48"],
    "6776": ["blaOXA-48"], "6804": ["blaVIM-1"], "6795": ["blaKPC-2"], "6813": ["blaKPC-2"],
    "6797": ["blaOXA-48"], "6838": ["blaOXA-48"], "6817": ["blaKPC-2"], "6855": ["blaKPC-2"],
    "6856": ["blaKPC-2"], "6888": ["blaKPC-3 (subtype only at 60x)", "blaKPC-2"],
    "6900": ["blaOXA-181"], "6901": ["blaKPC-2"], "6951": ["blaOXA-244"], "6970": ["blaKPC-2"],
    "6989": ["blaKPC-2"], "7000": ["blaNDM-1 (subtype not detected at 10x)"], "7067": ["blaVIM-1"],
    "7077": ["blaOXA-48"], "7084": ["blaKPC-2"], "7088": ["blaOXA-48"], "7089": ["blaOXA-48"],
    "7090": ["blaNDM-5"], "7083": ["blaOXA-48"], "7136": ["blaOXA-181"], "7176": ["blaOXA-48"],
    "7181": ["blaOXA-48", "blaVIM-1"], "7213": ["blaNDM-1"], "7221": ["blaNDM-5"],
    "7243": ["blaVIM-1"], "7285": ["blaKPC-2"], "7273": ["blaOXA-244"], "7316": ["blaOXA-181"],
    "7310": ["blaOXA-244"], "7370": ["blaNDM-5"], "7365": ["blaKPC-2"], "7418": ["blaKPC-2"],
    "7419": ["blaKPC-2"], "7387": ["blaOXA-244"], "7454": ["blaKPC-2"], "7476": ["blaKPC-2"],
    "7475": ["blaKPC-2"], "7488": ["blaOXA-48"], "7521": ["blaOXA-48"], "7503": ["blaNDM-5"]
}
S1B = {
    "5371": [], "5425": [], "5607": [], "5668": ["blaOXA-917, blaGES-5"], "5861": [], "6036": [],
    "6668": [], "6778": [], "7101": ["Imi-1"], "7320": []
}

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--catalog", required=True)
ap.add_argument("--vs-gold", default="bart_results_amr_vs_gold.tsv")
ap.add_argument("--sets", nargs="+", default=["carblot_pos", "carblot_neg"])
ap.add_argument("--depths", nargs="+", default=["10x", "20x"])
args = ap.parse_args()

sub = {}
with open(args.catalog, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r.get("allele"):
            sub[r["allele"]] = (r.get("subclass") or "").upper()

key = lambda s: re.sub(r"_\d+x$", "", s)
num = lambda s: (re.match(r"^isolate_(\d+)", s) or re.match(r"^(\d+)", s) or [None, None])[1]

def s1_genes(gl):
    """normalise the S1 'Carbapenemase gene (WGS)' entries to AMRFinderPlus symbols"""
    out = set()
    for g in gl:
        for t in re.split(r"[,;]", re.sub(r"\(.*?\)", "", g)):
            t = t.strip()
            if not t:
                continue
            t = re.sub(r"^(?i:bla)?", "", t)
            sym = "bla" + t[0].upper() + t[1:] if t else ""
            sym = re.sub(r"^bla(Imi|IMI)", "blaIMI", sym)
            if sym in sub:
                if "CARBAPENEM" in sub[sym]:
                    out.add(sym)
            elif re.match(r"^bla(KPC|NDM|VIM|IMP|OXA|GES|IMI)", sym):
                out.add(sym)            # not in catalogue: keep carbapenemase-like names (blaDHA etc. dropped)
    return out

# isolates present in the analysis (read folders)
seen = collections.defaultdict(set)        # isolate -> {set}
for s in args.sets:
    for d in args.depths:
        for p in sorted(Path(s, d, "bart").glob("*/")):
            seen[key(p.name)].add(s)

# reference: source and carbapenemase alleles
ref_src, ref_carb = {}, collections.defaultdict(set)
if Path(args.vs_gold).exists():
    with open(args.vs_gold, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ref_src[r["isolate"]] = r.get("gold_source", "")
            for a in r["gold_alleles"].split(";"):
                if a and ("CARBAPENEM" in sub.get(a, "") or r["gene_family"].endswith("_carbapenemase")):
                    ref_carb[r["isolate"]].add(a)
else:
    print(f"WARNING: {args.vs_gold} not found; reference columns empty", file=sys.stderr)

by_id = collections.defaultdict(list)
for iso in seen:
    by_id[num(iso)].append(iso)

rows, chosen = [], {}
for iso in sorted(seen, key=lambda s: (int(num(s) or 0), s)):
    i = num(iso)
    grp = "S1a" if i in S1A else "S1b" if i in S1B else "not_in_S1"
    s1 = s1_genes(S1A.get(i) or S1B.get(i) or [])
    ref = ref_carb.get(iso, set())
    dup = len(by_id[i]) > 1
    if grp != "not_in_S1":
        cand = sorted(by_id[i], key=lambda s: (s not in ref_src, s))   # prefer a run with reference
        chosen[i] = cand[0]
    match = "NA" if iso not in ref_src or grp == "not_in_S1" else ("yes" if ref == s1 else "no")
    rows.append(dict(isolate=iso, id=i, folder=";".join(sorted(seen[iso])), s1_group=grp,
                     reference=ref_src.get(iso, "none"), ref_carbapenemases=";".join(sorted(ref)) or "-",
                     s1_carbapenemases=";".join(sorted(s1)) or "-", genes_match_s1=match,
                     duplicate_runs="yes" if dup else "no",
                     used="yes" if chosen.get(i) == iso else "no"))

with open("isolate_classification.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)
with open("isolate_groups.tsv", "w") as fh:
    fh.write("isolate\tgroup\n")
    for r in rows:
        if r["used"] == "yes":
            fh.write(f'{r["isolate"]}\t{r["s1_group"]}\n')

# ---- report ----
P = lambda *a: print(*a)
c = collections.Counter((r["folder"], r["s1_group"]) for r in rows)
P("== analysed isolates by sequencing folder and S1 group")
for (f, g), n in sorted(c.items()):
    P(f"  {f:25s} {g:10s} {n}")
used = [r for r in rows if r["used"] == "yes"]
P(f"\n== used for the comparison: S1a {sum(r['s1_group']=='S1a' for r in used)}/{len(S1A)}, "
  f"S1b {sum(r['s1_group']=='S1b' for r in used)}/{len(S1B)}")
P("== S1 isolates not found in the analysis:", sorted(set(S1A) - set(by_id), key=int), sorted(set(S1B) - set(by_id), key=int))
P("== S1 isolates without reference assembly:", [r["isolate"] for r in used if r["reference"] == "none"])
P("== analysed isolates not in S1 (excluded from the comparison):")
for r in rows:
    if r["s1_group"] == "not_in_S1":
        P(f'  {r["isolate"]:45s} {r["folder"]:12s} ref carbapenemases: {r["ref_carbapenemases"]}')
P("== isolates with several runs (one used):")
for r in rows:
    if r["duplicate_runs"] == "yes":
        P(f'  {r["isolate"]:45s} used={r["used"]} reference={r["reference"]}')
P("== carbapenemase genes: reference assembly vs S1 'Carbapenemase gene (WGS)'")
for r in used:
    if r["genes_match_s1"] == "no":
        P(f'  {r["isolate"]:45s} {r["s1_group"]}  reference: {r["ref_carbapenemases"]:30s} S1: {r["s1_carbapenemases"]}')
P("\nwritten: isolate_classification.tsv, isolate_groups.tsv")
