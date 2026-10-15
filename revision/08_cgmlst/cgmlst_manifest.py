#!/usr/bin/env python3
"""Build the input sets for depth-dependent cgMLST with chewBBACA.

Takes the Flye/Medaka assemblies at every depth from the mlst tables (columns coverage, isolate,
assembly_path), keeps the isolates of Supplementary Table S1 (one run per isolate, from
isolate_groups.tsv), assigns each isolate to a cgMLST scheme from its species in Table S1 and links
the assemblies into <work>/<scheme>/inputs/<id>_<depth>.fasta (one chewBBACA run per scheme).

Writes <work>/manifest.tsv, <work>/training.tsv (assembly used for the Prodigal training file per
scheme: highest depth, lowest isolate number) and <work>/not_typed.tsv (species without a scheme).

Usage:
  python3 cgmlst_manifest.py --mlst-tables all_coverages_mlst_combined.tsv \
      --groups isolate_groups.tsv --work cgmlst_work
"""
import argparse, csv, re, sys, os, collections
from pathlib import Path

# species per isolate, Supplementary Table S1 (Pathogenwatch)
S1_SPECIES = {
    "3989": "Citrobacter freundii", "4704": "Citrobacter freundii", "4838": "Klebsiella pneumoniae",
    "4920": "Klebsiella pneumoniae", "5043": "Proteus mirabilis", "5244": "Klebsiella pneumoniae",
    "5305": "Providencia stuartii", "5371": "Enterobacter hormaechei", "5399": "Escherichia coli",
    "5425": "Citrobacter freundii", "5590": "Enterobacter hormaechei", "5607": "Enterobacter ludwigii",
    "5668": "Citrobacter freundii", "5860": "Klebsiella pneumoniae", "5861": "Klebsiella aerogenes",
    "5871": "Klebsiella variicola", "5969": "Citrobacter freundii", "5970": "Citrobacter freundii",
    "5976": "Enterobacter hormaechei", "6009": "Escherichia coli", "6014": "K. pneumoniae",
    "6035": "Klebsiella pneumoniae", "6036": "Enterobacter hormaechei", "6056": "Proteus mirabilis",
    "6057": "Escherichia coli", "6067": "Klebsiella michiganensis", "6068": "Escherichia coli",
    "6233": "Escherichia coli", "6243": "Klebsiella variicola", "6244": "Klebsiella pneumoniae",
    "6247": "Escherichia coli", "6287": "Klebsiella quasipneumoniae", "6312": "Escherichia coli",
    "6365": "Klebsiella variicola", "6378": "Klebsiella pneumoniae", "6394": "Klebsiella pneumoniae",
    "6403": "Escherichia coli", "6433": "Serratia marcescens", "6454": "Klebsiella variicola",
    "6455": "Klebsiella pneumoniae", "6466": "Escherichia coli", "6472": "Citrobacter freundii",
    "6482": "Enterobacter hormaechei", "6500": "Klebsiella pneumoniae",
    "6520": "Klebsiella pneumoniae", "6567": "Raoultella planticola", "6571": "Citrobacter freundii",
    "6601": "Citrobacter freundii", "6623": "Enterobacter hormaechei", "6668": "Klebsiella aerogenes",
    "6669": "Klebsiella pneumoniae", "6670": "Citrobacter freundii",
    "6701": "Klebsiella michiganensis", "6711": "Citrobacter europaeus", "6712": "Escherichia coli",
    "6726": "Klebsiella pneumoniae", "6735": "Citrobacter freundii", "6737": "Klebsiella pneumoniae",
    "6738": "Citrobacter freundii", "6748": "Escherichia coli", "6774": "Klebsiella pneumoniae",
    "6775": "Klebsiella aerogenes", "6776": "Citrobacter freundii", "6778": "Klebsiella aerogenes",
    "6779": "Klebsiella aerogenes", "6795": "Citrobacter amalonaticus", "6797": "Citrobacter freundii",
    "6804": "Enterobacter hormaechei", "6813": "Hafnia alvei", "6817": "Klebsiella pneumoniae",
    "6838": "Escherichia coli", "6855": "Klebsiella michiganensis", "6856": "Enterobacter hormaechei",
    "6888": "Klebsiella pneumoniae", "6900": "Citrobacter portucalensis",
    "6901": "Citrobacter freundii", "6951": "Escherichia coli", "6970": "Citrobacter freundii",
    "6989": "Klebsiella pneumoniae", "7000": "Proteus mirabilis", "7067": "Escherichia coli",
    "7077": "Serratia marcescens", "7083": "Escherichia coli", "7084": "E. roggenkampi",
    "7088": "Citrobacter koseri", "7089": "Klebsiella pneumoniae", "7090": "Citrobacter koseri",
    "7101": "Enterobacter bugandensis", "7136": "Klebsiella pneumoniae",
    "7176": "Klebsiella pneumoniae", "7181": "Klebsiella pneumoniae", "7213": "Proteus mirabilis",
    "7221": "Escherichia coli", "7243": "Enterobacter hormaechei", "7273": "Escherichia coli",
    "7285": "Citrobacter braakii", "7310": "Escherichia coli", "7316": "Klebsiella pneumoniae",
    "7320": "Klebsiella pneumoniae", "7365": "Klebsiella aerogenes", "7370": "Shigella dysenteriae",
    "7387": "Escherichia coli", "7418": "Citrobacter braakii", "7419": "Enterobacter hormaechei",
    "7454": "Escherichia coli", "7475": "Klebsiella pneumoniae", "7476": "Klebsiella pneumoniae",
    "7488": "Escherichia coli", "7503": "C. freundii", "7521": "Klebsiella pneumoniae"
}

# species -> scheme key (see schemes.tsv)
SCHEME_OF = {
    "Klebsiella pneumoniae": "kpsc", "K. pneumoniae": "kpsc", "Klebsiella variicola": "kpsc",
    "Klebsiella quasipneumoniae": "kpsc",
    "Escherichia coli": "ecoli", "Shigella dysenteriae": "ecoli",
    "Citrobacter freundii": "cfreundii", "C. freundii": "cfreundii", "Citrobacter braakii": "cfreundii",
    "Citrobacter europaeus": "cfreundii", "Citrobacter portucalensis": "cfreundii",
    "Enterobacter hormaechei": "ehorm",
    "Klebsiella michiganensis": "koxy",
    "Proteus mirabilis": "pmir",
    "Serratia marcescens": "smar",
}

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--mlst-tables", nargs="+", default=[])
ap.add_argument("--scan", action="store_true",
                help="also take assemblies directly from <base-dir>/polished_fasta_<depth>x/<isolate>_<depth>x.fasta")
ap.add_argument("--groups", required=True, help="isolate_groups.tsv from classify_isolates.py")
ap.add_argument("--work", default="cgmlst_work")
ap.add_argument("--base-dirs", nargs="+", default=["carblot_pos", "carblot_neg", "."],
                help="folders that relative assembly_path entries are resolved against, in this order")
args = ap.parse_args()

key = lambda s: re.sub(r"_\d+x$", "", s)
num = lambda s: (re.match(r"^isolate_(\d+)", s) or re.match(r"^(\d+)", s) or [None, None])[1]
dnum = lambda c: int(re.sub(r"\D", "", c) or 0)

groups = {}
with open(args.groups, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        groups[r["isolate"]] = r["group"]

paths = collections.defaultdict(set)          # (isolate, depth) -> paths
for t in args.mlst_tables:
    with open(t, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            iso, p = key(r["isolate"].strip()), (r.get("assembly_path") or "").strip()
            if iso in groups and p:
                paths[(iso, dnum(r["coverage"]))].add(p)

if args.scan:
    for b in args.base_dirs:
        for f in Path(b).glob("polished_fasta_*x/*.fasta"):
            m = re.fullmatch(r"(.+)_(\d+)x", f.stem)
            if m and m.group(1) in groups:
                paths[(m.group(1), int(m.group(2)))].add(str(f.resolve()))

work = Path(args.work)
rows, missing, conflicts, not_typed = [], [], [], set()
for (iso, d), ps in sorted(paths.items(), key=lambda x: (int(num(x[0][0])), -x[0][1])):
    i = num(iso)
    species = S1_SPECIES.get(i, "")
    scheme = SCHEME_OF.get(species)
    if not scheme:
        not_typed.add((i, iso, species or "unknown"))
        continue
    existing = sorted({str(c.resolve()) for p in ps for c in ([Path(p)] if Path(p).is_absolute() else
                       [Path(b) / p for b in args.base_dirs]) if c.is_file()})
    if not existing:
        missing.append(f"{iso} {d}x: {sorted(ps)}"); continue
    if len(existing) > 1:
        conflicts.append(f"{iso} {d}x: {existing}")
    src = Path(existing[0]).resolve()
    link = work / scheme / "inputs" / f"{i}_{d}x.fasta"
    link.parent.mkdir(parents=True, exist_ok=True)
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(src)
    rows.append(dict(id=i, isolate=iso, s1_group=groups[iso], species=species, scheme=scheme,
                     depth=d, genome=f"{i}_{d}x", path=str(src)))

print(f"isolate/depth combinations in the mlst tables matching isolate_groups.tsv: {len(paths)}; "
      f"assemblies found: {len(rows)}; listed but not found: {len(missing)}", file=sys.stderr)
if missing:
    print(f"not found, e.g. {missing[:3]}", file=sys.stderr)
if not rows:
    sys.exit("ERROR: no assemblies found; check the mlst tables, isolate_groups.tsv and --base-dirs")
with open(work / "manifest.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
with open(work / "not_typed.tsv", "w") as fh:
    fh.write("id\tisolate\tspecies\n" + "".join(f"{a}\t{b}\t{c}\n" for a, b, c in sorted(not_typed)))
train = {}
for r in sorted(rows, key=lambda r: (-r["depth"], int(r["id"]))):
    train.setdefault(r["scheme"], r)
with open(work / "training.tsv", "w") as fh:
    fh.write("scheme\tgenome\tpath\n" + "".join(f"{k}\t{v['genome']}\t{v['path']}\n" for k, v in sorted(train.items())))

c = collections.Counter((r["scheme"], r["depth"]) for r in rows)
iso_per = collections.defaultdict(set)
for r in rows:
    iso_per[r["scheme"]].add(r["id"])
depths = sorted({r["depth"] for r in rows}, reverse=True)
print("scheme      isolates  " + "  ".join(f"{d:>4}x" for d in depths))
for s in sorted(iso_per):
    print(f"{s:10s}  {len(iso_per[s]):8d}  " + "  ".join(f"{c.get((s, d), 0):5d}" for d in depths))
print(f"\nnot typed (no scheme): {len(not_typed)} isolates -> {work}/not_typed.tsv")
if missing:
    print(f"WARNING: {len(missing)} assemblies listed but not found, e.g. {missing[:3]}", file=sys.stderr)
if conflicts:
    print(f"WARNING: {len(conflicts)} isolate/depth with several assembly files (first used), e.g. {conflicts[:2]}", file=sys.stderr)
print(f"written: {work}/manifest.tsv, training.tsv, inputs per scheme")
