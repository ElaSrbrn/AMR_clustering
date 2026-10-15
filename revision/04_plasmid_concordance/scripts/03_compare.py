#!/usr/bin/env python3
"""Step 03: structural comparison of each test assembly with its reference plasmid.

For every row of work/pairs.tsv the whole test assembly is aligned to the reference plasmid
(minimap2 -cx asm5, primary alignments). The 'primary' contig is the test contig carrying the
carbapenemase gene; if several carry it, the one with the most reference sequence aligned is used.

Outputs: work/compare/<plasmid_id>/<assembler>_<depth>.paf
         results/stability_long.tsv   one row per plasmid x assembler x depth, with metrics and class
         results/breakpoints.tsv      every breakpoint on the primary contig, in reference coordinates
Usage  : $PY scripts/03_compare.py   (after `source config.sh`; needs minimap2 on PATH)

Classes (first matching rule wins):
  gene_absent            no contig carries the carbapenemase gene
  on_chromosome_contig   the gene sits on a contig >= CHROM_MIN (plasmid merged with or misplaced into the chromosome)
  concordant             |length diff| <= LEN_TOL, aligned fraction >= MIN_AF, no breakpoint, identity >= MIN_ID, circular
  concordant_not_circular  as concordant, but the test contig is not circular
  rearranged             >= 1 breakpoint that is a strand switch or a jump in reference coordinates
  expanded_duplication   reference fully covered, but part of it covered twice (duplication or overlapping contig ends)
  chimeric               reference fully covered, but > CHIMERA_UNALIGNED of the contig does not align to it
  internal_deletion      internal gap: reference sequence missing between two aligned blocks
  internal_insertion     internal gap: extra test sequence between two collinear blocks
  fragmented             reference coverage below MIN_AF on the primary contig, but other contigs cover the rest
  truncated              reference coverage below MIN_AF and not recovered by other contigs
  length_mismatch        structurally collinear but length differs by more than LEN_TOL (small indels)
  low_identity           everything fits except identity < MIN_ID
"""
import os
import subprocess
import sys

sys.path.insert(0, "lib")
from common import read_tsv, write_tsv, read_fasta, read_paf, union_len, circ_delta  # noqa: E402

E = os.environ
LEN_TOL = float(E.get("LEN_TOL", 0.001))
MIN_AF = float(E.get("MIN_AF", 0.99))
MIN_ID = float(E.get("MIN_ID", 0.999))
MIN_BLOCK = int(E.get("MIN_BLOCK", 1000))
BREAK_GAP = int(E.get("BREAK_GAP", 1000))
CHROM_MIN = int(E.get("CHROM_MIN", 1000000))
CHIMERA = float(E.get("CHIMERA_UNALIGNED", 0.05))
FRAG_MIN = int(E.get("FRAGMENT_MIN_BP", 2000))
THREADS = E.get("SLURM_CPUS_PER_TASK", "4")

refs = {r["plasmid_id"]: r for r in read_tsv("results/reference_plasmids.tsv")}
pairs = read_tsv("work/pairs.tsv")


def drop_contained(blocks):
    """Remove blocks whose query interval lies inside a longer block of the same contig."""
    keep = []
    for b in sorted(blocks, key=lambda x: -(x["qe"] - x["qs"])):
        if not any(k["qs"] <= b["qs"] and b["qe"] <= k["qe"] for k in keep):
            keep.append(b)
    return sorted(keep, key=lambda x: x["qs"])


def breakpoints(blocks, L, circular):
    """Walk blocks in query order and classify each discontinuity."""
    out = []
    for b1, b2 in zip(blocks, blocks[1:]):
        qgap = b2["qs"] - b1["qe"]
        if b1["strand"] != b2["strand"]:
            out.append(("rearrangement", b1, b2, None))
            continue
        if b1["strand"] == "+":
            tgap = circ_delta(b2["ts"] - b1["te"], L, circular)
        else:
            tgap = circ_delta(b1["ts"] - b2["te"], L, circular)
        if abs(tgap - qgap) <= BREAK_GAP:
            continue
        if abs(tgap) <= BREAK_GAP and qgap > BREAK_GAP:
            kind = "insertion"
        elif abs(qgap) <= BREAK_GAP and tgap > BREAK_GAP:
            kind = "deletion"
        else:
            kind = "rearrangement"
        out.append((kind, b1, b2, tgap - qgap))
    return out


long_rows, bp_rows = [], []
for p in pairs:
    pid, asm, dep = p["plasmid_id"], p["assembler"], p["depth"]
    ref = refs[pid]
    L = int(ref["length"])
    ref_circ = ref["circular"] == "true"
    row = {"plasmid_id": pid, "id": p["id"], "genes": p["genes"], "ref_length": L, "ref_rep_type": ref["rep_type"],
           "assembler": asm, "depth": dep}
    if p["status"] != "OK":
        row["class"] = "not_annotated"
        long_rows.append(row)
        continue
    paf = f"work/compare/{pid}/{asm}_{dep}.paf"
    if not os.path.exists(paf) or os.path.getsize(paf) == 0:
        with open(paf, "w") as fh:
            subprocess.run(["minimap2", "-cx", "asm5", "--secondary=no", "-t", THREADS,
                            f"work/compare/{pid}/ref.fasta", p["test_fasta"]],
                           stdout=fh, stderr=subprocess.DEVNULL, check=False)
    blocks = read_paf(paf, MIN_BLOCK, split_indel=BREAK_GAP)
    gene_contigs = [c for c in p["gene_contigs"].split(",") if c]
    qlens = {}
    circ = {}
    cpath = f"work/annot/{p['id']}/{asm}_{dep}/contigs.tsv"
    for c in read_tsv(cpath):
        qlens[c["contig"]] = int(c["length"])
        circ[c["contig"]] = c["circular"]
    by_q = {}
    for b in blocks:
        by_q.setdefault(b["q"], []).append(b)
    # whole-assembly coverage of the reference (any contig except chromosome-sized ones)
    cov_all = union_len([(b["ts"], b["te"]) for b in blocks if qlens.get(b["q"], 0) < CHROM_MIN])
    row["ref_af_all_nonchrom"] = round(cov_all / L, 4)
    row["n_gene_contigs"] = len(gene_contigs)
    if not gene_contigs:
        row.update({"class": "gene_absent", "ref_af_primary": 0})
        long_rows.append(row)
        continue
    prim = max(gene_contigs, key=lambda c: union_len([(b["ts"], b["te"]) for b in by_q.get(c, [])]))
    pb = drop_contained(by_q.get(prim, []))
    qlen = qlens.get(prim, 0)
    af = union_len([(b["ts"], b["te"]) for b in pb]) / L
    sum_ref = sum(b["te"] - b["ts"] for b in pb)
    dup_bp = sum_ref - union_len([(b["ts"], b["te"]) for b in pb])
    unaligned_q = 1 - union_len([(b["qs"], b["qe"]) for b in pb]) / qlen if qlen else 1
    ident = sum(b["nmatch"] for b in pb) / sum(b["alnlen"] for b in pb) if pb else 0
    bps = breakpoints(pb, L, ref_circ)
    n_rearr = sum(k == "rearrangement" for k, *_ in bps)
    n_ins = sum(k == "insertion" for k, *_ in bps)
    n_del = sum(k == "deletion" for k, *_ in bps)
    others = [c for c in by_q if c != prim and qlens.get(c, 0) < CHROM_MIN
              and union_len([(b["ts"], b["te"]) for b in by_q[c]]) >= FRAG_MIN]
    af_with_frag = union_len([(b["ts"], b["te"]) for c in [prim] + others for b in by_q.get(c, [])]) / L
    len_diff = (qlen - L) / L
    pcirc = circ.get(prim, "NA")
    rep_test = p["gene_contig_rep_types"].split("|")[gene_contigs.index(prim)] if p["gene_contig_rep_types"] else ""
    rs, ts_ = set(filter(None, ref["rep_type"].split(","))), set(filter(None, rep_test.split(",")))
    row.update({"primary_contig": prim, "primary_length": qlen, "primary_circular": pcirc,
                "length_diff_frac": round(len_diff, 5), "ref_af_primary": round(af, 4),
                "identity": round(ident, 5), "n_blocks": len(pb), "n_rearrangements": n_rearr,
                "n_insertions": n_ins, "n_deletions": n_del, "dup_bp": dup_bp,
                "unaligned_frac_primary": round(unaligned_q, 4), "n_fragment_contigs": len(others),
                "ref_af_with_fragments": round(af_with_frag, 4), "test_rep_type": rep_test,
                "replicon_match": ("NA" if not rs and not ts_ else ("yes" if rs == ts_ else ("partial" if rs & ts_ else "no")))})
    if qlen >= CHROM_MIN:
        cls = "on_chromosome_contig"
    elif abs(len_diff) <= LEN_TOL and af >= MIN_AF and not bps and ident >= MIN_ID:
        cls = "concordant" if pcirc == "true" else "concordant_not_circular"
    elif n_rearr:
        cls = "rearranged"
    elif af >= MIN_AF and dup_bp > BREAK_GAP:
        cls = "expanded_duplication"
    elif af >= MIN_AF and unaligned_q > CHIMERA:
        cls = "chimeric"
    elif n_del:
        cls = "internal_deletion"
    elif n_ins:
        cls = "internal_insertion"
    elif af < MIN_AF:
        cls = "fragmented" if af_with_frag >= MIN_AF or af_with_frag - af > 0.05 else "truncated"
    elif abs(len_diff) > LEN_TOL:
        cls = "length_mismatch"
    elif ident < MIN_ID:
        cls = "low_identity"
    else:
        cls = "other"
    row["class"] = cls
    long_rows.append(row)
    for kind, b1, b2, size in bps:
        for side, pos in (("left", b1["te"] if b1["strand"] == "+" else b1["ts"]),
                          ("right", b2["ts"] if b2["strand"] == "+" else b2["te"])):
            bp_rows.append({"plasmid_id": pid, "assembler": asm, "depth": dep, "class": cls, "primary_contig": prim,
                            "kind": kind, "side": side, "ref_pos": pos, "size_bp": "" if size is None else size,
                            "query_pos": b1["qe"] if side == "left" else b2["qs"]})

fields = ["plasmid_id", "id", "genes", "ref_length", "ref_rep_type", "assembler", "depth", "class",
          "primary_contig", "primary_length", "primary_circular", "length_diff_frac", "ref_af_primary",
          "identity", "n_blocks", "n_rearrangements", "n_insertions", "n_deletions", "dup_bp",
          "unaligned_frac_primary", "n_gene_contigs", "n_fragment_contigs", "ref_af_with_fragments",
          "ref_af_all_nonchrom", "test_rep_type", "replicon_match"]
write_tsv("results/stability_long.tsv", long_rows, fields)
write_tsv("results/breakpoints.tsv", bp_rows,
          ["plasmid_id", "assembler", "depth", "class", "primary_contig", "kind", "side", "ref_pos", "size_bp", "query_pos"])
print(f"Compared {len(long_rows)} plasmid x assembly pairs; {len(bp_rows) // 2} breakpoints.")
