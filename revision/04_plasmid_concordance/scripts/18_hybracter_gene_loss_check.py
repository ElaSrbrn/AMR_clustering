#!/usr/bin/env python3
"""Step 18: for Hybracter assemblies that lack the carbapenemase gene (class gene_absent in results/stability_long.tsv),
check where the plasmid was lost inside the Hybracter run (Supplementary Table S3c).

For each such assembly, the reference plasmid (work/compare/<plasmid>/ref.fasta) and its carbapenemase window
(work/context/<plasmid>/window_*.fasta from step 16) are aligned with Minimap2 (asm10) to
  1. Hybracter's intermediate Flye assembly   processing/assemblies/<run>/assembly.fasta
  2. the Plassembler output                   processing/plassembler/<run>/plassembler_plasmids.fasta
  3. the final assembly                       FINAL_OUTPUT/{complete,incomplete}/<run>_final.fasta
Reported: share of the query covered by alignment blocks >= 1 kb (union) and the longest single block.
Output: results/hybracter_gene_loss_check.tsv
Usage : HYB=/path/to/hybracter_depth_out_medaka python3 scripts/18_hybracter_gene_loss_check.py
Needs minimap2 on PATH (the mapping environment) and step 16 to have run.
"""
import csv
import glob
import os
import subprocess

H = os.environ.get("HYB", "/ictstr01/groups/hpc/urban_lab/projects/plasmid_project/hybracter_depth_out_medaka")


def cov(query, target):
    if not os.path.exists(target) or os.path.getsize(target) == 0:
        return "", "", "file missing"
    out = subprocess.run(["minimap2", "-c", "-x", "asm10", target, query], capture_output=True, text=True).stdout
    iv, ql, best = [], None, 0
    for line in out.splitlines():
        f = line.split("\t")
        ql, qs, qe = int(f[1]), int(f[2]), int(f[3])
        if qe - qs >= 1000:
            iv.append((qs, qe))
            best = max(best, qe - qs)
    if not iv:
        return "0", "0", ""
    iv.sort()
    tot, (cs, ce) = 0, iv[0]
    for s, e in iv[1:]:
        if s > ce:
            tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    tot += ce - cs
    return f"{100 * tot / ql:.0f}", f"{100 * best / ql:.0f}", ""


with open("results/stability_long.tsv") as fh:
    rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["assembler"] == "hybracter" and r["class"] == "gene_absent"]
out = []
for r in rows:
    pid, depth = r["plasmid_id"], r["depth"]
    iso = pid.rsplit("_", 1)[0]
    run = f"{iso}_{depth}"
    final = (glob.glob(f"{H}/{run}/FINAL_OUTPUT/*/{run}_final.fasta") or [""])[0]
    route = final.split("/FINAL_OUTPUT/")[1].split("/")[0] if final else "unknown"
    ref = f"work/compare/{pid}/ref.fasta"
    win = (sorted(glob.glob(f"work/context/{pid}/window_*.fasta")) or [""])[0]
    for label, path in [("intermediate Flye assembly", f"{H}/{run}/processing/assemblies/{run}/assembly.fasta"),
                        ("Plassembler output", f"{H}/{run}/processing/plassembler/{run}/plassembler_plasmids.fasta"),
                        ("final assembly", final)]:
        pc, pb, note = cov(ref, path) if path else ("", "", "file missing")
        wc, wb, _ = cov(win, path) if (path and win) else ("", "", "")
        out.append({"plasmid_id": pid, "depth": depth, "hybracter_route": route, "file": label,
                    "plasmid_covered_pct": pc, "plasmid_best_block_pct": pb,
                    "window_covered_pct": wc, "window_best_block_pct": wb, "note": note})
        print("\t".join(str(v) for v in out[-1].values()))
with open("results/hybracter_gene_loss_check.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, ["plasmid_id", "depth", "hybracter_route", "file", "plasmid_covered_pct",
                            "plasmid_best_block_pct", "window_covered_pct", "window_best_block_pct", "note"], delimiter="\t")
    w.writeheader()
    w.writerows(out)
print(f"written: results/hybracter_gene_loss_check.tsv ({len(rows)} Hybracter assemblies without the carbapenemase gene)")
