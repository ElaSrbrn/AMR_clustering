#!/usr/bin/env python3
"""Step 07: IS elements and repeats on the selected reference plasmids, and where the breakpoints fall.

Inputs : work/indepth/<plasmid_id>/{isfinder_ref,isfinder_test,self}.blast6, results/breakpoints.tsv
Outputs: work/indepth/<plasmid_id>/is_ref.bed, repeats.bed, is_test.tsv
         results/is_elements_reference.tsv     IS elements on each selected reference plasmid
         results/is_elements_test_contigs.tsv  IS elements on the divergent test contigs
         results/breakpoints_annotated.tsv     each breakpoint with distance to the nearest IS element / repeat
         results/breakpoint_is_summary.tsv     per assembler x depth: breakpoints within IS_WINDOW of an IS element or repeat
"""
import os
import sys
from collections import defaultdict

sys.path.insert(0, "lib")
from common import read_tsv, write_tsv, union_len  # noqa: E402

E = os.environ
IS_MIN_ID = float(E.get("IS_MIN_ID", 90))
IS_MIN_COV = float(E.get("IS_MIN_COV", 80))
REP_LEN = int(E.get("REPEAT_MIN_LEN", 500))
REP_ID = float(E.get("REPEAT_MIN_ID", 90))
WINDOW = int(E.get("IS_WINDOW", 1000))
F6 = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"]


def blast(path):
    out = []
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 12:
                r = dict(zip(F6, p))
                for k in ("pident", "bitscore"):
                    r[k] = float(r[k])
                for k in ("length", "qstart", "qend", "sstart", "send", "qlen", "slen"):
                    r[k] = int(r[k])
                out.append(r)
    return out


def is_hits(rows):
    """Filter ISfinder hits and keep the best non-overlapping hit per region (greedy by bitscore), per query."""
    keep = defaultdict(list)
    good = [r for r in rows if r["pident"] >= IS_MIN_ID and 100 * r["length"] / r["slen"] >= IS_MIN_COV]
    for r in sorted(good, key=lambda x: -x["bitscore"]):
        s, e = sorted((r["qstart"], r["qend"]))
        if any(min(e, ke) - max(s, ks) > 0.5 * (e - s) for ks, ke, _ in keep[r["qseqid"]]):
            continue
        keep[r["qseqid"]].append((s, e, r))
    return keep


def background(intervals, L):
    """Fraction of the plasmid lying within WINDOW of an IS element or repeat (expected rate for random breakpoints)."""
    iv = []
    for s, e in intervals:
        a, b = s - WINDOW, e + WINDOW
        if a < 0:
            iv += [(0, b), (L + a, L)]
        elif b > L:
            iv += [(a, L), (0, b - L)]
        else:
            iv.append((a, b))
    return min(1.0, union_len([(max(0, a), min(L, b)) for a, b in iv]) / L) if L else 0


def dist(pos, intervals, L):
    best = None
    for s, e in intervals:
        if s <= pos <= e:
            return 0
        d = min(abs(pos - s), abs(pos - e))
        d = min(d, L - d) if L else d
        best = d if best is None else min(best, d)
    return best


refs = {r["plasmid_id"]: r for r in read_tsv("results/reference_plasmids.tsv")}
sel = [l.strip() for l in open("work/indepth_list.txt") if l.strip()] if os.path.exists("work/indepth_list.txt") else []
is_ref_rows, is_test_rows, feat = [], [], {}
for pid in sel:
    d = f"work/indepth/{pid}"
    L = int(refs[pid]["length"])
    ish = is_hits(blast(f"{d}/isfinder_ref.blast6")).get(pid, [])
    with open(f"{d}/is_ref.bed", "w") as fb:
        for s, e, r in sorted(ish):
            strand = "+" if r["sstart"] < r["send"] else "-"
            fb.write(f"{pid}\t{s - 1}\t{e}\t{r['sseqid']}\t{r['pident']:.1f}\t{strand}\n")
            is_ref_rows.append({"plasmid_id": pid, "start": s, "end": e, "is_element": r["sseqid"], "identity": r["pident"],
                                "coverage": round(100 * r["length"] / r["slen"], 1), "strand": strand})
    reps = []
    for r in blast(f"{d}/self.blast6"):
        if r["qstart"] == r["sstart"] and r["qend"] == r["send"]:
            continue  # diagonal
        if r["length"] >= REP_LEN and r["pident"] >= REP_ID:
            reps.append(tuple(sorted((r["qstart"], r["qend"]))))
    reps = sorted(set(reps))
    with open(f"{d}/repeats.bed", "w") as fb:
        for s, e in reps:
            fb.write(f"{pid}\t{s - 1}\t{e}\trepeat\n")
    feat[pid] = ([(s, e) for s, e, _ in ish], reps, L)
    for q, hits in is_hits(blast(f"{d}/isfinder_test.blast6")).items():
        for s, e, r in sorted(hits):
            is_test_rows.append({"plasmid_id": pid, "test_contig": q, "start": s, "end": e, "is_element": r["sseqid"],
                                 "identity": r["pident"], "coverage": round(100 * r["length"] / r["slen"], 1)})
    write_tsv(f"{d}/is_test.tsv", [x for x in is_test_rows if x["plasmid_id"] == pid],
              ["plasmid_id", "test_contig", "start", "end", "is_element", "identity", "coverage"])

write_tsv("results/is_elements_reference.tsv", is_ref_rows,
          ["plasmid_id", "start", "end", "is_element", "identity", "coverage", "strand"])
write_tsv("results/is_elements_test_contigs.tsv", is_test_rows,
          ["plasmid_id", "test_contig", "start", "end", "is_element", "identity", "coverage"])

ann = []
summ = defaultdict(lambda: [0, 0, 0, 0.0])
bg = {pid: background(f[0] + f[1], f[2]) for pid, f in feat.items()}
for b in read_tsv("results/breakpoints.tsv"):
    if b["plasmid_id"] not in feat:
        continue
    isi, reps, L = feat[b["plasmid_id"]]
    pos = int(b["ref_pos"])
    di, dr = dist(pos, isi, L), dist(pos, reps, L)
    near_is = di is not None and di <= WINDOW
    near_any = near_is or (dr is not None and dr <= WINDOW)
    ann.append({**b, "dist_is": "" if di is None else di, "dist_repeat": "" if dr is None else dr,
                "near_is": "yes" if near_is else "no", "near_is_or_repeat": "yes" if near_any else "no"})
    k = (b["assembler"], b["depth"])
    summ[k][0] += 1
    summ[k][1] += near_is
    summ[k][2] += near_any
    summ[k][3] += bg[b["plasmid_id"]]
write_tsv("results/breakpoints_annotated.tsv", ann,
          ["plasmid_id", "assembler", "depth", "class", "primary_contig", "kind", "side", "ref_pos", "size_bp",
           "query_pos", "dist_is", "dist_repeat", "near_is", "near_is_or_repeat"])
srows = [{"assembler": a, "depth": d, "breakpoint_ends": n, "near_is": i, "near_is_or_repeat": r,
          "prop_near_is_or_repeat": round(r / n, 3) if n else "",
          "expected_prop_random": round(x / n, 3) if n else ""} for (a, d), (n, i, r, x) in sorted(summ.items())]
write_tsv("results/breakpoint_is_summary.tsv", srows,
          ["assembler", "depth", "breakpoint_ends", "near_is", "near_is_or_repeat", "prop_near_is_or_repeat",
           "expected_prop_random"])
print(f"IS elements on {len(sel)} reference plasmids: {len(is_ref_rows)}; on divergent test contigs: {len(is_test_rows)}")
tot = sum(v[0] for v in summ.values())
near = sum(v[2] for v in summ.values())
exp = sum(v[3] for v in summ.values())
if tot:
    print(f"Breakpoint ends within {WINDOW} bp of an IS element or repeat: {near}/{tot} ({near / tot:.1%}); "
          f"expected for random positions on the same plasmids: {exp / tot:.1%}")
    print("Breakpoint ends of one breakpoint are not independent; test enrichment per breakpoint or per plasmid, not per end.")
