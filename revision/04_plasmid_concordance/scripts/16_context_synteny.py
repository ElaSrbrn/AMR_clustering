#!/usr/bin/env python3
"""Step 16: synteny of the carbapenemase genetic context, reference vs test assemblies.
For every carbapenemase gene on each reference plasmid (work/compare/<plasmid>/ref.fasta): locate it with AMRFinderPlus,
extract gene +/- FLANK bp from the circular reference, align the window with minimap2 to the whole Flye and Hybracter
assembly of the isolate at each depth (paths from work/manifest.tsv), and classify: syntenic (one alignment covers
>= MIN_COV of the window), syntenic_across_contig_start (only after joining contig end to start; valid only for
circular contigs), not_syntenic.
Outputs: results/context_synteny_long.tsv, results/context_synteny_summary.tsv, results/context_synteny_files.tsv
"""
import argparse, csv, glob, os, subprocess, sys, tempfile
from collections import defaultdict


def read_fasta(path):
    recs, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def write_fasta(path, recs):
    with open(path, "w") as fh:
        for n, s in recs:
            fh.write(f">{n}\n")
            for i in range(0, len(s), 80):
                fh.write(s[i:i + 80] + "\n")


def circ_window(seq, start, end, flank):
    L = len(seq)
    glen = end - start + 1
    if glen + 2 * flank >= L:
        s0 = (start - 1 - (L - glen) // 2) % L
        return (seq + seq)[s0:s0 + L], s0, L, True
    s0 = (start - 1 - flank) % L
    n = glen + 2 * flank
    return (seq + seq)[s0:s0 + n], s0, n, False


def carbapenemase_genes(ref, out, db):
    if not (os.path.exists(out) and os.path.getsize(out) > 0):
        cmd = ["amrfinder", "-n", ref, "-o", out, "--threads", "2"]
        if db:
            cmd += ["-d", db]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    genes = []
    with open(out) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if "CARBAPENEM" not in (r.get("Subclass") or "").upper():
                continue
            sym = r.get("Element symbol") or r.get("Gene symbol") or "NA"
            genes.append((sym, r["Contig id"], int(r["Start"]), int(r["Stop"]), r["Strand"]))
    return genes


def best_hit(window_fa, recs, tmp, double, preset):
    tfa = os.path.join(tmp, "target.fa")
    write_fasta(tfa, [(n, s + s if double else s) for n, s in recs])
    r = subprocess.run(["minimap2", "-c", "-x", preset, "--secondary=no", tfa, window_fa],
                       capture_output=True, text=True, check=True)
    best = None
    for line in r.stdout.splitlines():
        f = line.split("\t")
        qlen, qs, qe = int(f[1]), int(f[2]), int(f[3])
        cov = (qe - qs) / qlen
        ident = int(f[9]) / int(f[10]) if int(f[10]) else 0.0
        if best is None or cov > best["cov"]:
            best = {"cov": cov, "ident": ident, "contig": f[5], "strand": f[4]}
    return best


def read_manifest(path):
    out, skipped = defaultdict(list), defaultdict(int)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["assembler"] not in ("flye", "hybracter"):
                continue
            if r.get("status", "OK") != "OK" or not os.path.exists(r["path"]):
                skipped[(r["assembler"], r["depth"])] += 1
                continue
            out[r["id"]].append((r["assembler"], r["depth"], r["path"]))
    return out, skipped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--flank", type=int, default=10000)
    ap.add_argument("--min-cov", type=float, default=0.95)
    ap.add_argument("--preset", default="asm10")
    ap.add_argument("--manifest", default="work/manifest.tsv")
    a = ap.parse_args()
    db = os.environ.get("AMRFINDER_DB", "")
    plasmid_dirs = sorted(d for d in glob.glob("work/compare/*") if os.path.exists(os.path.join(d, "ref.fasta")))
    print(f"reference plasmids found: {len(plasmid_dirs)} (expected 93)")
    os.makedirs("results", exist_ok=True)
    long_rows, file_rows, no_gene, no_asm = [], [], [], []
    man, skipped = read_manifest(a.manifest)
    print(f"manifest: {sum(len(v) for v in man.values())} Flye/Hybracter assemblies for {len(man)} isolates; "
          f"skipped (status not OK or file missing): {dict(skipped)}")
    for d in plasmid_dirs:
        pid = os.path.basename(d)
        ref = os.path.join(d, "ref.fasta")
        refseq = dict(read_fasta(ref))
        cdir = os.path.join("work/context", pid)
        os.makedirs(cdir, exist_ok=True)
        genes = carbapenemase_genes(ref, os.path.join(cdir, "amrfinder.tsv"), db)
        if not genes:
            no_gene.append(pid)
            continue
        iso = pid.rsplit("_", 1)[0]
        tests = sorted(man.get(iso, []))
        if not tests:
            no_asm.append(pid)
        for asm, depth, f in tests:
            file_rows.append({"plasmid_id": pid, "isolate": iso, "assembler": asm, "depth": depth, "path": f})
        with tempfile.TemporaryDirectory() as tmp:
            for gi, (sym, ctg, st, en, strand) in enumerate(genes, 1):
                win, s0, wlen, whole = circ_window(refseq[ctg], st, en, a.flank)
                wfa = os.path.join(tmp, "window.fa")
                write_fasta(wfa, [(f"{pid}_{sym}", win)])
                write_fasta(os.path.join(cdir, f"window_{gi}_{sym}.fasta"), [(f"{pid}_{sym}_{s0}", win)])
                for asm, depth, f in tests:
                    recs = read_fasta(f)
                    h = best_hit(wfa, recs, tmp, False, a.preset)
                    status = "syntenic" if h and h["cov"] >= a.min_cov else None
                    h2 = None
                    if status is None:
                        h2 = best_hit(wfa, recs, tmp, True, a.preset)
                        status = "syntenic_across_contig_start" if h2 and h2["cov"] >= a.min_cov else "not_syntenic"
                    use = h2 if status == "syntenic_across_contig_start" else h
                    long_rows.append({
                        "plasmid_id": pid, "gene": sym, "gene_start": st, "gene_end": en, "gene_strand": strand,
                        "window_bp": wlen, "window_is_whole_plasmid": "yes" if whole else "no",
                        "assembler": asm, "depth": depth, "status": status,
                        "best_cov_pct": f"{100 * use['cov']:.1f}" if use else "0.0",
                        "best_identity_pct": f"{100 * use['ident']:.2f}" if use else "",
                        "best_contig": use["contig"] if use else ""})

    def wtsv(path, rows):
        if not rows:
            open(path, "w").close()
            return
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    wtsv("results/context_synteny_long.tsv", long_rows)
    wtsv("results/context_synteny_files.tsv", file_rows)
    summ = defaultdict(lambda: defaultdict(int))
    for r in long_rows:
        k = (r["assembler"], r["depth"])
        summ[k]["genes_assessed"] += 1
        summ[k][r["status"]] += 1
    srows = [{"assembler": asm, "depth": depth, "genes_assessed": c["genes_assessed"], "syntenic": c["syntenic"],
              "syntenic_across_contig_start": c["syntenic_across_contig_start"], "not_syntenic": c["not_syntenic"]}
             for (asm, depth), c in sorted(summ.items())]
    wtsv("results/context_synteny_summary.tsv", srows)
    print(f"flank {a.flank} bp, min coverage {a.min_cov}, minimap2 preset {a.preset}")
    print(f"plasmids without a carbapenemase gene found by AMRFinderPlus: {len(no_gene)} {' '.join(no_gene)}")
    print(f"plasmids without any Flye/Hybracter assembly in the manifest: {len(no_asm)} {' '.join(no_asm)}")
    per = defaultdict(int)
    for r in file_rows:
        per[r["plasmid_id"]] += 1
    print(f"test assemblies used: {len(file_rows)}; plasmids with fewer than 8 (2 assemblers x 4 depths): "
          f"{sum(1 for v in per.values() if v < 8)}")
    for r in srows:
        print("\t".join(str(v) for v in r.values()))


if __name__ == "__main__":
    sys.exit(main())
