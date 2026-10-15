"""Shared helpers for the plasmid contig stability pipeline (Python standard library only)."""
import csv
import os
import re

COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def read_fasta(path):
    """Return list of (name, header, sequence). name = first word of the header."""
    recs, name, header, seq = [], None, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, header, "".join(seq)))
                header = line[1:]
                name = header.split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        recs.append((name, header, "".join(seq)))
    return recs


def write_fasta(path, recs, width=80):
    with open(path, "w") as f:
        for name, seq in recs:
            f.write(f">{name}\n")
            for i in range(0, len(seq), width):
                f.write(seq[i:i + width] + "\n")


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def read_tsv(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows, fields):
    with open(path, "w") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def col(row, *names, default=""):
    for n in names:
        if n in row and row[n] is not None:
            return row[n]
    return default


# ---------------------------------------------------------------- circularity
_CIRC_RE = re.compile(r"circular=(\w+)", re.I)


def circular_from_header(header):
    m = _CIRC_RE.search(header)
    if not m:
        return None
    return m.group(1).lower() in ("true", "yes", "y", "1")


def circular_from_side_tables(fasta_path):
    """Look for Flye assembly_info.txt or Hybracter per_contig_stats next to the FASTA. Returns {contig: bool}."""
    d = os.path.dirname(os.path.abspath(fasta_path))
    out = {}
    cand = [os.path.join(d, "assembly_info.txt")]
    # Medaka-polished Flye: .../polished_flye_<depth>/<name>_medaka/consensus.fasta
    #   -> Flye output  .../flye_<depth>/<name>_flye/assembly_info.txt
    parent, leaf = os.path.split(d)
    if leaf.endswith("_medaka") and os.path.basename(parent).startswith("polished_flye_"):
        flye_dir = os.path.join(os.path.dirname(parent), os.path.basename(parent)[len("polished_"):],
                                leaf[:-len("_medaka")] + "_flye")
        cand.append(os.path.join(flye_dir, "assembly_info.txt"))
    cand += [os.path.join(d, x) for x in os.listdir(d) if x.endswith("per_contig_stats.tsv")]
    for p in cand:
        if not os.path.exists(p):
            continue
        with open(p) as f:
            lines = [l.rstrip("\n").split("\t") for l in f if l.strip()]
        if not lines:
            continue
        hdr = [h.lstrip("#").strip().lower() for h in lines[0]]
        ci = next((i for i, h in enumerate(hdr) if "circ" in h), None)
        if ci is None:
            continue
        ni = 0
        for row in lines[1:]:
            if len(row) > ci:
                out[row[ni]] = row[ci].strip().lower() in ("y", "yes", "true", "1", "+")
    return out


def hybracter_incomplete_circularity(fasta_path, recs, max_rel_diff=0.005):
    """Hybracter 'incomplete' output has no circularity flag. Its contigs are the (polished, renamed)
    contigs of the internal Flye assembly, whose assembly_info.txt records circularity:
    <sample>/FINAL_OUTPUT/incomplete/<sample>_final.fasta -> <sample>/processing/assemblies/<sample>/assembly_info.txt
    Contigs are matched by length (closest Flye contig within max_rel_diff, each Flye contig used once).
    Returns {contig: bool} for matched contigs only."""
    d = os.path.dirname(os.path.abspath(fasta_path))
    if os.path.basename(d) != "incomplete" or os.path.basename(os.path.dirname(d)) != "FINAL_OUTPUT":
        return {}
    sample_dir = os.path.dirname(os.path.dirname(d))
    sample = os.path.basename(sample_dir)
    info = os.path.join(sample_dir, "processing", "assemblies", sample, "assembly_info.txt")
    if not os.path.exists(info):
        return {}
    flye = []
    with open(info) as f:
        hdr = None
        for line in f:
            p = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = [h.lstrip("#").strip().lower() for h in p]
                li = next(i for i, h in enumerate(hdr) if h.startswith("length"))
                ci = next(i for i, h in enumerate(hdr) if "circ" in h)
                continue
            if len(p) > max(li, ci):
                flye.append((int(p[li]), p[ci].strip().upper() in ("Y", "YES", "TRUE", "+")))
    used, out = set(), {}
    for name, _, seq in sorted(recs, key=lambda r: -len(r[2])):
        L = len(seq)
        best = None
        for j, (fl, circ) in enumerate(flye):
            if j in used:
                continue
            rd = abs(fl - L) / max(L, 1)
            if rd <= max_rel_diff and (best is None or rd < best[0]):
                best = (rd, j, circ)
        if best:
            used.add(best[1])
            out[name] = best[2]
    return out


def contig_table(fasta_path):
    """List of dicts: contig, length, circular ('true'/'false'/'NA')."""
    recs = read_fasta(fasta_path)
    side = None
    if not any(circular_from_header(h) is not None for _, h, _ in recs):
        side = circular_from_side_tables(fasta_path)
        side.update(hybracter_incomplete_circularity(fasta_path, recs))
    rows = []
    for name, header, seq in recs:
        c = circular_from_header(header)
        if c is None:
            if side is None:
                side = circular_from_side_tables(fasta_path)
            c = side.get(name)
        rows.append({"contig": name, "length": len(seq),
                     "circular": "NA" if c is None else ("true" if c else "false")})
    return rows


# ---------------------------------------------------------------- annotation parsers
def amr_hits(path):
    """AMRFinderPlus table -> list of dicts with contig, gene, subclass, start, stop, strand."""
    out = []
    for r in read_tsv(path):
        out.append({
            "contig": col(r, "Contig id", "Contig ID", "Contig"),
            "gene": col(r, "Element symbol", "Gene symbol"),
            "subclass": col(r, "Subclass", "Element subclass"),
            "start": col(r, "Start"),
            "stop": col(r, "Stop"),
            "strand": col(r, "Strand"),
            "name": col(r, "Element name", "Sequence name"),
            "method": col(r, "Method"),
        })
    return out


# Carbapenemase families whose partial calls at a contig end (gene split across the start/end of a circular
# contig) AMRFinderPlus reports only at family level, without the CARBAPENEM subclass.
_CARB_PREFIX = ("blaKPC", "blaNDM", "blaVIM", "blaIMP", "blaGES", "blaIMI", "blaNMC", "blaSME", "blaFRI", "blaGIM", "blaSPM")


def carb_family(hit):
    """Family key used to match carbapenemase calls between assemblies (allele-independent)."""
    name, sym = hit.get("name", ""), hit["gene"]
    if "OXA-48 family" in name or "OXA-48-like" in name:
        return "blaOXA-48-like"
    for pre in _CARB_PREFIX:
        if sym.startswith(pre):
            return pre
    return sym


def is_carbapenemase(hit):
    if "CARBAPENEM" in hit["subclass"].upper():
        return True
    # partial call at a contig end of a known carbapenemase family (e.g. blaOXA-232 spanning the origin of a small circular plasmid)
    partial_end = "PARTIAL_CONTIG_END" in hit.get("method", "").upper() or "PARTIAL" in hit.get("method", "").upper()
    return partial_end and (carb_family(hit) in _CARB_PREFIX or carb_family(hit) == "blaOXA-48-like")


def mob_types(path):
    """MOB-typer --multi table -> {contig: {rep_type, relaxase, mobility}}."""
    out = {}
    for r in read_tsv(path):
        sid = col(r, "sample_id").split()[0]
        sid = sid.split(":")[-1]
        out[sid] = {"rep_type": col(r, "rep_type(s)"),
                    "relaxase": col(r, "relaxase_type(s)"),
                    "mobility": col(r, "predicted_mobility")}
    return out


# ---------------------------------------------------------------- PAF
_CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def _split_block(b, cigar, big):
    """Split one PAF block at insertions/deletions >= big (from the cg:Z CIGAR) into collinear sub-blocks."""
    ops = [(int(n), o) for n, o in _CIG.findall(cigar)]
    m_total = sum(n for n, o in ops if o in "M=X") or 1
    rate = b["nmatch"] / m_total
    fwd = b["strand"] == "+"
    q = b["qs"]
    t = b["ts"] if fwd else b["te"]          # reference walks backwards on the minus strand
    out, cur = [], None
    for n, o in ops:
        if cur is None:
            cur = {"qs": q, "tpos": t, "m": 0, "aln": 0}
        if o in "M=X":
            q += n
            t += n if fwd else -n
            cur["m"] += n
            cur["aln"] += n
        elif o == "I":
            if n >= big:
                out.append((cur, q, t)); cur = None
            else:
                cur["aln"] += n
            q += n
        elif o in "DN":
            if n >= big:
                out.append((cur, q, t)); cur = None
            else:
                cur["aln"] += n
            t += n if fwd else -n
    if cur is not None:
        out.append((cur, q, t))
    blocks = []
    for c, qe, tend in out:
        if c["m"] == 0:
            continue
        ts, te = (c["tpos"], tend) if fwd else (tend, c["tpos"])
        blocks.append({**b, "qs": c["qs"], "qe": qe, "ts": ts, "te": te,
                       "nmatch": int(round(c["m"] * rate)), "alnlen": c["aln"]})
    return blocks


def read_paf(path, min_block=0, split_indel=None):
    """Primary alignments only. Each block: q, qlen, qs, qe, strand, t, tlen, ts, te, nmatch, alnlen.
    If split_indel is set and the PAF has cg:Z tags (minimap2 -c), blocks are split at indels >= split_indel."""
    blocks = []
    if not os.path.exists(path):
        return blocks
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            tags = p[12:]
            if any(t.startswith("tp:A:") and t != "tp:A:P" for t in tags):
                continue
            b = {"q": p[0], "qlen": int(p[1]), "qs": int(p[2]), "qe": int(p[3]), "strand": p[4],
                 "t": p[5], "tlen": int(p[6]), "ts": int(p[7]), "te": int(p[8]),
                 "nmatch": int(p[9]), "alnlen": int(p[10])}
            cg = next((t[5:] for t in tags if t.startswith("cg:Z:")), None)
            parts = _split_block(b, cg, split_indel) if (split_indel and cg) else [b]
            blocks += [x for x in parts if x["qe"] - x["qs"] >= min_block]
    return blocks


def union_len(intervals):
    tot, cur_s, cur_e = 0, None, None
    for s, e in sorted(intervals):
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                tot += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        tot += cur_e - cur_s
    return tot


def circ_delta(d, L, circular):
    """Smallest absolute offset allowing wrap-around on a circular reference."""
    if not circular or L <= 0:
        return d
    return min((d, d - L, d + L), key=abs)
