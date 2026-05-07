#!/usr/bin/env python3
import argparse
import csv
import re
import shutil
from pathlib import Path

CARB_PATTERNS = [
    re.compile(r'\bKPC(?:-\d+)?\b', re.I),
    re.compile(r'\bNDM(?:-\d+)?\b', re.I),
    re.compile(r'\bVIM(?:-\d+)?\b', re.I),
    re.compile(r'\bIMP(?:-\d+)?\b', re.I),
    re.compile(r'\bIMI(?:-\d+)?\b', re.I),
    re.compile(r'\bSME(?:-\d+)?\b', re.I),
    re.compile(r'\bGIM(?:-\d+)?\b', re.I),
    re.compile(r'\bSIM(?:-\d+)?\b', re.I),
    re.compile(r'\bDIM(?:-\d+)?\b', re.I),
    re.compile(r'\bSPM(?:-\d+)?\b', re.I),
    re.compile(r'\bAIM(?:-\d+)?\b', re.I),
    re.compile(r'\bNMC-A\b', re.I),
    re.compile(r'\bOXA-48(?:-LIKE)?\b', re.I),
    re.compile(r'\bOXA-181\b', re.I),
    re.compile(r'\bOXA-232\b', re.I),
    re.compile(r'\bOXA-244\b', re.I),
    re.compile(r'\bOXA-162\b', re.I),
    re.compile(r'\bOXA-163\b', re.I),
    re.compile(r'\bOXA-370\b', re.I),
    re.compile(r'\bOXA-23(?:-LIKE)?\b', re.I),
    re.compile(r'\bOXA-24(?:/40|-LIKE)?\b', re.I),
    re.compile(r'\bOXA-58(?:-LIKE)?\b', re.I),
]

FAMILY_ORDER = ["KPC", "NDM", "VIM", "IMP", "IMI", "SME", "GIM", "SIM", "DIM", "SPM", "AIM", "NMC-A", "OXA"]

def norm(x):
    return (x or "").strip()

def pick_col(fieldnames, candidates):
    if not fieldnames:
        return None
    lowered = {f.lower(): f for f in fieldnames}
    for c in candidates:
        if c.lower() in lowered:
            return lowered[c.lower()]
    for c in candidates:
        for f in fieldnames:
            if c.lower() in f.lower():
                return f
    return None

def infer_family(text):
    u = text.upper().replace("BLA", "")
    for fam in FAMILY_ORDER:
        if fam in u:
            return fam
    return ""

def is_carbapenemase(symbol, name, subclass, subtype):
    text = " ".join([symbol, name, subclass, subtype]).upper().replace("BLA", "")
    if any(p.search(text) for p in CARB_PATTERNS):
        return True
    if "CARBAPENEM" in text:
        return True
    return False

def fasta_iter(path: Path):
    header = None
    seq = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
 yield header, "".join(seq)
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if header is not None:
            yield header, "".join(seq)

def write_single_fasta(outpath: Path, header: str, seq: str):
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with outpath.open("w") as out:
        out.write(f">{header}\n")
        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")

def parse_contig_report(path: Path):
    contigs = {}
    if not path.exists():
        return contigs
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        f = reader.fieldnames or []
        contig_col = pick_col(f, ["contig_id", "contig", "id"])
        mol_col = pick_col(f, ["molecule_type", "type"])
        cluster_col = pick_col(f, ["primary_cluster_id", "cluster_id"])
        secondary_col = pick_col(f, ["secondary_cluster_id"])
        rep_col = pick_col(f, ["rep_type(s)", "rep_type", "replicon_type"])
        mobility_col = pick_col(f, ["predicted_mobility", "mobility"])
        nn_col = pick_col(f, ["mash_neighbor_identification", "mash_nearest_neighbor"])
        size_col = pick_col(f, ["size"])
        circ_col = pick_col(f, ["circularity_status"])
        for row in reader:
            contig_id = norm(row.get(contig_col, "")) if contig_col else ""
            if not contig_id:
                continue
            molecule_type = norm(row.get(mol_col, "")) if mol_col else ""
            annotation = molecule_type.lower() if molecule_type else "unknown"
            if annotation not in {"plasmid", "chromosome"}:
                annotation = "plasmid" if norm(row.get(cluster_col, "")) else "unknown"
            contigs[contig_id] = {
                "annotation": annotation,
                "primary_cluster_id": norm(row.get(cluster_col, "")) if cluster_col else "",
                "secondary_cluster_id": norm(row.get(secondary_col, "")) if secondary_col else "",
                "rep_type": norm(row.get(rep_col, "")) if rep_col else "",
                "predicted_mobility": norm(row.get(mobility_col, "")) if mobility_col else "",
                "mash_neighbor": norm(row.get(nn_col, "")) if nn_col else "",
                "size": norm(row.get(size_col, "")) if size_col else "",
                "circularity_status": norm(row.get(circ_col, "")) if circ_col else "",
            }
    return contigs

def parse_amrfinder(path: Path):
    hits = []
    if not path.exists() or path.stat().st_size == 0:
        return hits
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        f = reader.fieldnames or []
        symbol_col = pick_col(f, ["Element symbol", "Gene symbol", "gene_symbol"])
        name_col = pick_col(f, ["Element name", "Gene name", "Protein name", "Sequence name"])
        contig_col = pick_col(f, ["Contig id", "contig_id", "contig"])
        subclass_col = pick_col(f, ["Subclass"])
        subtype_col = pick_col(f, ["Subtype"])
        scope_col = pick_col(f, ["Scope"])
        type_col = pick_col(f, ["Type"])
        method_col = pick_col(f, ["Method"])
        start_col = pick_col(f, ["Start"])
        stop_col = pick_col(f, ["Stop"])
        strand_col = pick_col(f, ["Strand"])
        cov_col = pick_col(f, ["% Coverage of reference"])
       ident_col = pick_col(f, ["% Identity to reference"])
        closest_acc_col = pick_col(f, ["Closest reference accession"])
        closest_name_col = pick_col(f, ["Closest reference name"])
        for row in reader:
            symbol = norm(row.get(symbol_col, "")) if symbol_col else ""
            name = norm(row.get(name_col, "")) if name_col else ""
            subclass = norm(row.get(subclass_col, "")) if subclass_col else ""
            subtype = norm(row.get(subtype_col, "")) if subtype_col else ""
            if not is_carbapenemase(symbol, name, subclass, subtype):
                continue
            hits.append({
                "gene_symbol": symbol,
                "gene_name": name,
                "contig_id": norm(row.get(contig_col, "")) if contig_col else "",
                "subclass": subclass,
                "subtype": subtype,
                "scope": norm(row.get(scope_col, "")) if scope_col else "",
                "type": norm(row.get(type_col, "")) if type_col else "",
                "method": norm(row.get(method_col, "")) if method_col else "",
                "start": norm(row.get(start_col, "")) if start_col else "",
                "stop": norm(row.get(stop_col, "")) if stop_col else "",
                "strand": norm(row.get(strand_col, "")) if strand_col else "",
                "pct_ref_coverage": norm(row.get(cov_col, "")) if cov_col else "",
                "pct_identity": norm(row.get(ident_col, "")) if ident_col else "",
                "closest_reference_accession": norm(row.get(closest_acc_col, "")) if closest_acc_col else "",
                "closest_reference_name": norm(row.get(closest_name_col, "")) if closest_name_col else "",
            })
    return hits

def find_assembly(assembly_root: Path, sample: str):
    for ext in (".fasta", ".fa", ".fna"):
        p = assembly_root / f"{sample}{ext}"
        if p.exists():
            return p
    matches = list(assembly_root.glob(f"{sample}.*"))
    for p in matches:
        if p.suffix.lower() in {".fasta", ".fa", ".fna"}:
           return None

def build_summary(annotation_root: Path, assembly_root: Path, coverage: str):
    summary_dir = annotation_root / "summaries"
    plasmid_dir = annotation_root / "carbapenemase_encoding_plasmids"
    summary_dir.mkdir(parents=True, exist_ok=True)
    plasmid_dir.mkdir(parents=True, exist_ok=True)

    tsv_path = summary_dir / f"carbapenemase_summary_{coverage}.tsv"
    csv_path = summary_dir / f"carbapenemase_summary_{coverage}.csv"

    fieldnames = [
        "sample", "coverage", "assembly_path", "amrfinder_file", "mob_contig_report",
        "gene_symbol", "gene_name", "carbapenemase_family", "contig_id",
        "start", "stop", "strand", "scope", "type", "subtype", "subclass", "method",
        "pct_ref_coverage", "pct_identity", "closest_reference_accession", "closest_reference_name",
        "contig_annotation", "primary_cluster_id", "secondary_cluster_id", "replicon_type",
        "predicted_mobility", "circularity_status", "mash_neighbor", "contig_size",
        "plasmid_exported", "plasmid_export_path"
    ]

    rows = []
    exported = set()

    for sample_dir in sorted(p for p in annotation_root.iterdir() if p.is_dir() and p.name not in {"summaries", "carbapenemase_encoding_plasmids", "selected_plasmids"}):
        sample = sample_dir.name
        amr = sample_dir / "amrfinder" / f"{sample}_assembly_amrfinder.tsv"
        mob = sample_dir / "mob_recon" / "contig_report.txt"
        assembly = find_assembly(assembly_root, sample)
        contig_info = parse_contig_report(mob) if mob.exists() else {}
        carb_hits = parse_amrfinder(amr)
        if not carb_hits:
            continue

        seqs = {}
        if assembly and assembly.exists():
            seqs = {h.split()[0]: s for h, s in fasta_iter(assembly)}

        for hit in carb_hits:
            cid = hit["contig_id"]
            cinfo = contig_info.get(cid, {})
            annotation = cinfo.get("annotation", "unknown")
            family = infer_family(" ".join([hit["gene_symbol"], hit["gene_name"], hit["subtype"], hit["subclass"]]))
            export_path = ""
            exported_flag = "no"

            if annotation == "plasmid" and cid in seqs:
                cluster = cinfo.get("primary_cluster_id") or "unclustered"
                family_dir = family or "other"
                outpath = plasmid_dir / family_dir / coverage / f"{sample}__{cid}__{cluster}.fasta"
                if outpath not in exported:
                    write_single_fasta(outpath, cid, seqs[cid])
                    exported.add(outpath)
                export_path = str(outpath)
                exported_flag = "yes"

            rows.append({
                "sample": sample,
                "coverage": coverage,
                "assembly_path": str(assembly) if assembly else "",
                "amrfinder_file": str(amr),
                "mob_contig_report": str(mob) if mob.exists() else "",
                "gene_symbol": hit["gene_symbol"],
                "gene_name": hit["gene_name"],
                "carbapenemase_family": family,
                "contig_id": cid,
                "start": hit["start"],
                "stop": hit["stop"],
                "strand": hit["strand"],
"scope": hit["scope"],
                "type": hit["type"],
                "subtype": hit["subtype"],
                "subclass": hit["subclass"],
                "method": hit["method"],
                "pct_ref_coverage": hit["pct_ref_coverage"],
                "pct_identity": hit["pct_identity"],
                "closest_reference_accession": hit["closest_reference_accession"],
                "closest_reference_name": hit["closest_reference_name"],
                "contig_annotation": annotation,
                "primary_cluster_id": cinfo.get("primary_cluster_id", ""),
                "secondary_cluster_id": cinfo.get("secondary_cluster_id", ""),
                "replicon_type": cinfo.get("rep_type", ""),
                "predicted_mobility": cinfo.get("predicted_mobility", ""),
                "circularity_status": cinfo.get("circularity_status", ""),
                "mash_neighbor": cinfo.get("mash_neighbor", ""),
                "contig_size": cinfo.get("size", ""),
                "plasmid_exported": exported_flag,
                "plasmid_export_path": export_path,
            })

    with tsv_path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    with csv_path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} carbapenemase rows to {tsv_path}")
    print(f"Wrote CSV to {csv_path}")
    print(f"Exported plasmid contigs to {plasmid_dir}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Build carbapenemase summary and export carbapenemase-encoding plasmid contigs.")
    ap.add_argument("--annotation-root", required=True, type=Path)
    ap.add_argument("--assembly-root", required=True, type=Path)
    ap.add_argument("--coverage", required=True)
    args = ap.parse_args()
    build_summary(args.annotation_root, args.assembly_root, args.coverage)
