#!/usr/bin/env python3
"""Step 02: define the reference carbapenemase plasmids and pair them with every test assembly.

Inputs : work/manifest.tsv, work/annot/<id>/<assembler>_<depth>/{contigs,amr,mob}.tsv
Outputs: results/reference_carbapenemase_loci.tsv  every carbapenemase gene in the references (plasmid or chromosome)
         results/reference_plasmids.tsv            the plasmids that enter the stability analysis
         work/compare/<plasmid_id>/ref.fasta        reference plasmid sequence
         work/pairs.tsv                             one row per reference plasmid x test assembly
Usage  : $PY scripts/02_collect.py   (from the pipeline folder, after `source config.sh`)
"""
import os
import sys

sys.path.insert(0, "lib")
from common import read_tsv, write_tsv, read_fasta, write_fasta, amr_hits, is_carbapenemase, carb_family, mob_types  # noqa: E402

CHROM_MIN = int(os.environ.get("CHROM_MIN", 1000000))

man = [r for r in read_tsv("work/manifest.tsv") if r["status"] == "OK"]
by_id = {}
for r in man:
    by_id.setdefault(r["id"], {})[(r["assembler"], r["depth"])] = r["path"]


def load(iso, asm, dep):
    d = f"work/annot/{iso}/{asm}_{dep}"
    if not os.path.exists(os.path.join(d, "DONE")):
        return None
    contigs = {c["contig"]: c for c in read_tsv(os.path.join(d, "contigs.tsv"))}
    amr = amr_hits(os.path.join(d, "amr.tsv"))
    mob = mob_types(os.path.join(d, "mob.tsv"))
    return contigs, amr, mob


loci, plasmids, pairs = [], [], []
missing_ref = []
for iso in sorted(by_id):
    ref = load(iso, "autocycler", "60x")
    if ref is None:
        missing_ref.append(iso)
        continue
    contigs, amr, mob = ref
    carb = [h for h in amr if is_carbapenemase(h)]
    by_contig = {}
    for h in carb:
        by_contig.setdefault(h["contig"], []).append(h)
        L = int(contigs.get(h["contig"], {}).get("length", 0))
        loci.append({"id": iso, "gene": h["gene"], "ref_contig": h["contig"], "start": h["start"], "stop": h["stop"],
                     "strand": h["strand"], "contig_length": L,
                     "circular": contigs.get(h["contig"], {}).get("circular", "NA"),
                     "location": "chromosome" if L >= CHROM_MIN else "plasmid",
                     "rep_type": mob.get(h["contig"], {}).get("rep_type", "")})
    ref_seqs = {n: s for n, _, s in read_fasta(by_id[iso][("autocycler", "60x")])}
    for c, hits in by_contig.items():
        L = int(contigs[c]["length"])
        if L >= CHROM_MIN:
            continue  # chromosomal genes are reported in loci but not part of the plasmid stability analysis
        pid = f"{iso}_{c}"
        genes = sorted({h["gene"] for h in hits})
        fams = {carb_family(h) for h in hits}
        m = mob.get(c, {})
        plasmids.append({"plasmid_id": pid, "id": iso, "ref_contig": c, "genes": ",".join(genes),
                         "n_gene_copies": len(hits), "length": L, "circular": contigs[c]["circular"],
                         "rep_type": m.get("rep_type", ""), "relaxase": m.get("relaxase", ""),
                         "mobility": m.get("mobility", "")})
        os.makedirs(f"work/compare/{pid}", exist_ok=True)
        write_fasta(f"work/compare/{pid}/ref.fasta", [(c, ref_seqs[c])])
        for (asm, dep), path in sorted(by_id[iso].items()):
            if asm == "autocycler":
                continue
            t = load(iso, asm, dep)
            if t is None:
                pairs.append({"plasmid_id": pid, "id": iso, "assembler": asm, "depth": dep, "test_fasta": path,
                              "genes": ",".join(genes), "gene_contigs": "", "status": "NOT_ANNOTATED"})
                continue
            tcontigs, tamr, tmob = t
            gc = sorted({h["contig"] for h in tamr if h["gene"] in genes or (is_carbapenemase(h) and carb_family(h) in fams)})
            pairs.append({"plasmid_id": pid, "id": iso, "assembler": asm, "depth": dep, "test_fasta": path,
                          "genes": ",".join(genes), "gene_contigs": ",".join(gc),
                          "gene_contig_lengths": ",".join(tcontigs.get(x, {}).get("length", "") for x in gc),
                          "gene_contig_rep_types": "|".join(tmob.get(x, {}).get("rep_type", "") for x in gc),
                          "status": "OK"})

os.makedirs("results", exist_ok=True)
write_tsv("results/reference_carbapenemase_loci.tsv", loci,
          ["id", "gene", "ref_contig", "start", "stop", "strand", "contig_length", "circular", "location", "rep_type"])
write_tsv("results/reference_plasmids.tsv", plasmids,
          ["plasmid_id", "id", "ref_contig", "genes", "n_gene_copies", "length", "circular", "rep_type", "relaxase", "mobility"])
write_tsv("work/pairs.tsv", pairs,
          ["plasmid_id", "id", "assembler", "depth", "test_fasta", "genes", "gene_contigs", "gene_contig_lengths",
           "gene_contig_rep_types", "status"])

n_chrom = sum(1 for l in loci if l["location"] == "chromosome")
print(f"Isolates with a reference: {len(by_id) - len(missing_ref)}  (missing/not annotated: {' '.join(missing_ref) or 'none'})")
print(f"Carbapenemase genes in references: {len(loci)}  (plasmid: {len(loci) - n_chrom}, chromosome: {n_chrom})")
print(f"Reference plasmids for the stability analysis: {len(plasmids)}")
print(f"Pairs (plasmid x test assembly): {len(pairs)}  not annotated: {sum(p['status'] != 'OK' for p in pairs)}")
iso_without = sorted(set(by_id) - set(missing_ref) - {l['id'] for l in loci})
if iso_without:
    print(f"WARNING: no carbapenemase gene found in the reference of: {' '.join(iso_without)}")
