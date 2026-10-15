#!/usr/bin/env python3
"""Build BART's AMR KMA index from a local AMRFinderPlus database (AMR_CDS[.fa]),
so read-level (BART) and assembly-level (AMRFinderPlus) AMR use the same database version.

BART expects template names "<gene_symbol> <description_with_underscores>" (it prints the
gene from the first word and the description with '_' replaced by spaces). bart-update -amr
takes the gene symbol from '|' field 6 and the description from field 8 of the NCBI header.

Usage:
  python3 build_bart_amr_index.py --show  <AMRDB_dir>             # inspect header fields first
  python3 build_bart_amr_index.py <AMRDB_dir> <out.fasta> [--symbol-field 6] [--desc-field 8]
Then:  kma index -i <out.fasta> -o <bart_db>/indexes/amr
Fields are 1-based.
"""
import sys, re, glob, os, collections

def find_cds(dbdir):
    for name in ("AMR_CDS.fa", "AMR_CDS.fasta", "AMR_CDS"):
        p = os.path.join(dbdir, name)
        if os.path.isfile(p):
            return p
    hits = glob.glob(os.path.join(dbdir, "AMR_CDS*"))
    if hits:
        return hits[0]
    sys.exit(f"ERROR: no AMR_CDS file in {dbdir}")

args = [a for a in sys.argv[1:]]
if not args:
    sys.exit(__doc__)

if args[0] == "--show":
    cds = find_cds(args[1])
    print(f"file: {cds}")
    n = 0
    with open(cds) as fh:
        for line in fh:
            if line.startswith(">"):
                print("\nheader:", line.rstrip()[:300])
                for i, f in enumerate(line[1:].rstrip().split("|"), 1):
                    print(f"  field {i}: {f[:100]}")
                n += 1
                if n == 3:
                    break
    sys.exit(0)

sym_f, desc_f = 6, 8
if "--symbol-field" in args:
    sym_f = int(args[args.index("--symbol-field") + 1])
if "--desc-field" in args:
    desc_f = int(args[args.index("--desc-field") + 1])
dbdir, out = args[0], args[1]
cds = find_cds(dbdir)

n, bad, examples, symbols = 0, 0, [], collections.Counter()
with open(cds) as fh, open(out, "w") as o:
    for line in fh:
        if line.startswith(">"):
            fields = line[1:].rstrip("\n").split("|")
            if len(fields) < max(sym_f, desc_f):
                bad += 1
                sym, desc = fields[0].split()[0], "unknown"
            else:
                sym = fields[sym_f - 1].strip()
                desc = fields[desc_f - 1].strip()
            sym = re.sub(r"\s+", "_", sym) or "unknown"
            desc = re.sub(r"\s+", "_", desc) or "unknown"
            o.write(f">{sym} {desc}\n")
            symbols[sym] += 1
            n += 1
            if len(examples) < 5:
                examples.append(f"{sym} {desc}")
        else:
            o.write(line)

print(f"source:    {cds}")
print(f"templates: {n}  (headers with too few fields: {bad})")
print("examples:  " + " | ".join(examples))
carb = sorted(s for s in symbols if re.match(r"bla(KPC|NDM|VIM|IMP|OXA-48)", s))
print(f"carbapenemase symbols found: {len(carb)} e.g. {carb[:6]}")
if bad:
    sys.exit(f"ERROR: {bad} headers have fewer than {max(sym_f, desc_f)} fields; check field numbers with --show")
if n == 0 or not carb:
    sys.exit("ERROR: no carbapenemase gene symbols found; check --symbol-field with --show")
print(f"next: kma index -i {out} -o <bart_db>/indexes/amr")
