#!/usr/bin/env python3
"""Step 08: one folder per selected plasmid for the collaborator's synteny figures (pyGenomeViz).

Reads selection/synteny.tsv (include == yes). For each plasmid writes synteny_packages/<plasmid_id>/:
  reference.fasta               Autocycler 60x reference plasmid
  reference.gbff / .gff3        Bakta annotation (if step 06 ran)
  <assembler>_<depth>.fasta     primary test contig, oriented to the reference strand; circular contigs rotated
                                to start at the position matching reference position 1
  is_elements.bed, repeats.bed  ISfinder hits and self-repeats on the reference (step 07)
  breakpoints.tsv               breakpoints of all test contigs in reference coordinates
  carbapenemase.bed             carbapenemase gene(s) on the reference
  README.txt                    what each file is, plus the stability class of every test contig
"""
import os
import shutil
import sys

sys.path.insert(0, "lib")
from common import read_tsv, write_tsv, read_fasta, write_fasta, revcomp, read_paf  # noqa: E402

if not os.path.exists("selection/synteny.tsv"):
    sys.exit("selection/synteny.tsv not found: copy results/proposal_synteny.tsv there and edit 'include'.")
sel = [r for r in read_tsv("selection/synteny.tsv") if r.get("include", "").strip().lower() == "yes"]
refs = {r["plasmid_id"]: r for r in read_tsv("results/reference_plasmids.tsv")}
pairs = {(p["plasmid_id"], p["assembler"], p["depth"]): p for p in read_tsv("work/pairs.tsv")}
long = read_tsv("results/stability_long.tsv")
bps = read_tsv("results/breakpoints.tsv")
loci = read_tsv("results/reference_carbapenemase_loci.tsv")
cache = {}

for s in sel:
    pid = s["plasmid_id"]
    if pid not in refs:
        print(f"skip {pid}: not a reference plasmid")
        continue
    out = f"synteny_packages/{pid}"
    os.makedirs(out, exist_ok=True)
    ref = refs[pid]
    rseq = read_fasta(f"work/compare/{pid}/ref.fasta")[0][2]
    write_fasta(f"{out}/reference.fasta", [(pid, rseq)])
    bk = f"work/indepth/{pid}/bakta"
    for ext in ("gbff", "gff3"):
        f = f"{bk}/{pid}.{ext}"
        if os.path.exists(f):
            shutil.copy(f, f"{out}/reference.{ext}")
    for name in ("is_ref.bed", "repeats.bed"):
        f = f"work/indepth/{pid}/{name}"
        if os.path.exists(f):
            shutil.copy(f, f"{out}/{'is_elements.bed' if name == 'is_ref.bed' else name}")
    with open(f"{out}/carbapenemase.bed", "w") as fb:
        for l in loci:
            if l["id"] == ref["id"] and l["ref_contig"] == ref["ref_contig"]:
                fb.write(f"{pid}\t{int(l['start']) - 1}\t{l['stop']}\t{l['gene']}\t0\t{l['strand']}\n")
    write_tsv(f"{out}/breakpoints.tsv", [b for b in bps if b["plasmid_id"] == pid],
              ["assembler", "depth", "class", "primary_contig", "kind", "side", "ref_pos", "size_bp", "query_pos"])
    readme = [f"Plasmid {pid}  genes: {ref['genes']}  length: {ref['length']} bp  circular: {ref['circular']}  "
              f"replicon: {ref['rep_type']}", f"Selected because: {s.get('reason', '')}", "",
              "Test contigs (assembler_depth: class, contig length, circular, aligned fraction of reference):"]
    for r in sorted((x for x in long if x["plasmid_id"] == pid), key=lambda x: (x["assembler"], x["depth"])):
        tag = f"{r['assembler']}_{r['depth']}"
        readme.append(f"  {tag}: {r['class']}  {r.get('primary_length', '')} bp  circular={r.get('primary_circular', '')}  "
                      f"aligned={r.get('ref_af_primary', '')}")
        prim = r.get("primary_contig")
        if not prim:
            continue
        fa = pairs[(pid, r["assembler"], r["depth"])]["test_fasta"]
        if fa not in cache:
            cache[fa] = {n: q for n, _, q in read_fasta(fa)}
        seq = cache[fa].get(prim)
        if not seq or len(seq) >= 1000000:
            continue
        blocks = [b for b in read_paf(f"work/compare/{pid}/{tag}.paf") if b["q"] == prim]
        if blocks:
            b = max(blocks, key=lambda x: x["qe"] - x["qs"])
            qlen = len(seq)
            if b["strand"] == "-":
                seq = revcomp(seq)
                qs = qlen - b["qe"]
            else:
                qs = b["qs"]
            if r.get("primary_circular") == "true":
                q0 = (qs - b["ts"]) % qlen
                seq = seq[q0:] + seq[:q0]
        write_fasta(f"{out}/{tag}.fasta", [(f"{tag}__{prim}", seq)])
    readme += ["", "Files:",
               "  reference.fasta / .gbff / .gff3   reference plasmid and Bakta annotation",
               "  <assembler>_<depth>.fasta         primary test contig (oriented; circular contigs rotated to the reference start)",
               "  is_elements.bed / repeats.bed     ISfinder hits and self-repeats on the reference (0-based BED)",
               "  carbapenemase.bed                 carbapenemase gene(s) on the reference",
               "  breakpoints.tsv                   breakpoints in reference coordinates (kind: insertion/deletion/rearrangement)",
               "Orientation and rotation use the longest alignment block; contigs with rearrangements can still show crossing links."]
    with open(f"{out}/README.txt", "w") as f:
        f.write("\n".join(readme) + "\n")
print(f"Packaged {len(sel)} plasmids in synteny_packages/  (tar -czf synteny_packages.tar.gz synteny_packages)")
