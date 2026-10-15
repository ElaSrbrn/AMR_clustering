#!/usr/bin/env python3
"""Step 09: FASTA sets for the Mash/DCJ (Pling) clustering, 60x reference vs a test assembler/depth.

Writes clustering_input/reference_60x/<plasmid_id>.fasta and clustering_input/<assembler>_<depth>/<plasmid_id>.fasta
(primary contig of each comparison; plasmids without a primary contig, or with the gene on a chromosome-sized
contig, are listed in clustering_input/<assembler>_<depth>/excluded.tsv instead).
Usage: $PY scripts/09_export_for_clustering.py [assembler_depth ...]     default: flye_40x
"""
import os
import sys

sys.path.insert(0, "lib")
from common import read_tsv, write_tsv, read_fasta, write_fasta  # noqa: E402

sets = sys.argv[1:] or ["flye_40x"]
CHROM_MIN = int(os.environ.get("CHROM_MIN", 1000000))
refs = read_tsv("results/reference_plasmids.tsv")
os.makedirs("clustering_input/reference_60x", exist_ok=True)
for r in refs:
    seq = read_fasta(f"work/compare/{r['plasmid_id']}/ref.fasta")[0][2]
    write_fasta(f"clustering_input/reference_60x/{r['plasmid_id']}.fasta", [(r["plasmid_id"], seq)])
pairs = {(p["plasmid_id"], p["assembler"], p["depth"]): p for p in read_tsv("work/pairs.tsv")}
long = read_tsv("results/stability_long.tsv")
cache = {}
for tag in sets:
    asm, dep = tag.split("_")
    out = f"clustering_input/{tag}"
    os.makedirs(out, exist_ok=True)
    excl, n = [], 0
    for r in long:
        if r["assembler"] != asm or r["depth"] != dep:
            continue
        prim = r.get("primary_contig")
        if not prim or r["class"] == "on_chromosome_contig":
            excl.append({"plasmid_id": r["plasmid_id"], "class": r["class"]})
            continue
        fa = pairs[(r["plasmid_id"], asm, dep)]["test_fasta"]
        if fa not in cache:
            cache[fa] = {k: s for k, _, s in read_fasta(fa)}
        write_fasta(f"{out}/{r['plasmid_id']}.fasta", [(r["plasmid_id"], cache[fa][prim])])
        n += 1
    write_tsv(f"{out}/excluded.tsv", excl, ["plasmid_id", "class"])
    print(f"{tag}: {n} plasmid contigs exported, {len(excl)} excluded (see {out}/excluded.tsv)")
print(f"reference_60x: {len(refs)} plasmids exported")
