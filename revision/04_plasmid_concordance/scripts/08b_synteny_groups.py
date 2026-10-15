#!/usr/bin/env python3
"""Step 08b: group packages for synteny comparisons BETWEEN reference plasmids (sets 2 and 3).
Reads synteny_groups.tsv and the per-plasmid packages of step 08 (synteny_packages/<plasmid_id>/).
Writes synteny_groups/<group>/: rotated .gbff and .fasta per plasmid (start 1 kb upstream of the carbapenemase gene,
gene on + strand), carbapenemase.tsv, pairwise_minimap2.paf (consecutive plasmids), README.txt.
"""
import os, shutil, subprocess, sys
from collections import defaultdict
sys.path.insert(0, "lib")
from common import read_tsv, write_fasta  # noqa: E402

FLANK = 1000
groups = defaultdict(list)
for r in read_tsv("synteny_groups.tsv"):
    groups[r["group"]].append(r)
try:
    from Bio import SeqIO
except ImportError:
    sys.exit("Biopython is needed (activate an env that has it, e.g. bakta_clean)")


def bed_first(path):
    for l in open(path):
        p = l.rstrip("\n").split("\t")
        if len(p) >= 6:
            return int(p[1]), int(p[2]), p[3], p[5]
    return None


mm2 = shutil.which("minimap2") or "/ictstr01/home/haicu/ela.sauerborn/mamba/envs/autocycler/bin/minimap2"
for g, rows in groups.items():
    rows.sort(key=lambda r: int(r["order"]))
    out = f"synteny_groups/{g}"
    os.makedirs(out, exist_ok=True)
    made, info = [], []
    for r in rows:
        pid = r["plasmid_id"]
        pk = f"synteny_packages/{pid}"
        gb, bed = f"{pk}/reference.gbff", f"{pk}/carbapenemase.bed"
        if not (os.path.exists(gb) and os.path.exists(bed)):
            print(f"[{g}] skip {pid}: {gb} or {bed} missing")
            continue
        rec = SeqIO.read(gb, "genbank")
        s, e, gene, strand = bed_first(bed)
        L = len(rec)
        if strand == "-":
            rec = rec.reverse_complement(id=rec.id, name=rec.name, description=rec.description, annotations=rec.annotations)
            s, e = L - e, L - s
        start = (s - FLANK) % L
        rec = rec[start:] + rec[:start]   # features spanning the cut are dropped by Biopython
        rec.annotations["molecule_type"] = "DNA"
        name = f"{r['order']}_{pid}"
        SeqIO.write(rec, f"{out}/{name}.gbff", "genbank")
        write_fasta(f"{out}/{name}.fasta", [(name, str(rec.seq))])
        made.append(name)
        info.append({"file": name, "label": r["label"], "gene": gene, "start": FLANK + 1, "end": FLANK + (e - s), "length": L})
    with open(f"{out}/carbapenemase.tsv", "w") as f:
        f.write("file\tlabel\tgene\tstart\tend\tplasmid_length\n")
        for i in info:
            f.write(f"{i['file']}\t{i['label']}\t{i['gene']}\t{i['start']}\t{i['end']}\t{i['length']}\n")
    if os.path.exists(mm2):
        with open(f"{out}/pairwise_minimap2.paf", "w") as f:
            for a, b in zip(made, made[1:]):
                subprocess.run([mm2, "-x", "asm5", "-c", f"{out}/{b}.fasta", f"{out}/{a}.fasta"],
                               stdout=f, stderr=subprocess.DEVNULL, check=False)
    else:
        print(f"[{g}] minimap2 not found: pairwise_minimap2.paf not written")
    readme = [f"Group: {g}", "", "Plasmids, top to bottom:"] + [
        f"  {i['file']}: {i['label']}; {i['gene']} at {i['start']}-{i['end']} of {i['length']} bp" for i in info] + [
        "", "All sequences are rotated to start 1 kb upstream of the carbapenemase gene and oriented with the gene on the + strand,",
        "so the gene lines up at the left of every track. Features spanning the original start position are dropped.",
        "Plot: one track per plasmid in this order, full length, genes coloured by category (carbapenemase red, other AMR orange,",
        "IS/transposase grey, replication/conjugation blue, other light grey), links from pairwise_minimap2.paf (or pgv-blast /",
        "pgv-mmseqs on the .gbff files) with an identity colour scale. A second version zoomed to the carbapenemase region",
        "(about 10 kb either side) shows the immediate genetic context.",
        "Please report the pyGenomeViz version used."]
    with open(f"{out}/README.txt", "w") as f:
        f.write("\n".join(readme) + "\n")
    print(f"[{g}] {len(made)} plasmids -> {out}/")
