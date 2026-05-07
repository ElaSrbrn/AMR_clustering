#!/usr/bin/env python3
import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

def sanitize(name: str) -> str:
    name = name.strip()
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
    return name or "unknown"

def split_replicons(rep_value: str) -> List[str]:
    if not rep_value:
        return []
    raw = re.split(r"[;,]", rep_value)
    groups = []
    for item in raw:
        item = item.strip()
        if item and item != "-":
            groups.append(item)
    return sorted(set(groups))

def read_fasta_length(path: Path) -> int:
    total = 0
    if path.suffix == ".gz":
        import gzip
        opener = gzip.open
    else:
        opener = open
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.startswith(">"):
                total += len(line.strip())
    return total

def choose_summary_file(annotation_dir: Path) -> Path:
    summaries = sorted((annotation_dir / "summaries").glob("carbapenemase_summary_*.tsv"))
    if not summaries:
        raise FileNotFoundError(f"no carbapenemase summary found under {annotation_dir}/summaries")
    return summaries[-1]

def link_or_copy(src: Path, dst: Path, mode: str) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    if mode == "copy":
        shutil.copy2(src, dst)
    else:
        rel = os.path.relpath(src.resolve(), start=dst.parent.resolve())
        os.symlink(rel, dst)

def run_checked(cmd: List[str], stdout_path: Path = None, stderr_path: Path = None) -> None:
    if stdout_path:
        stdout_handle = open(stdout_path, "w")
    else:
        stdout_handle = subprocess.DEVNULL
    if stderr_path:
        stderr_handle = open(stderr_path, "w")
    else:
        stderr_handle = subprocess.DEVNULL
    try:
        proc = subprocess.run(cmd, stdout=stdout_handle, stderr=stderr_handle, check=False, text=True)
        if proc.returncode != 0:
            raise RuntimeError("command failed: " + " ".join(cmd))
    finally:
        if stdout_path:
            stdout_handle.close()

      if stderr_path:
            stderr_handle.close()

def mash_sketch(mash_cmd: str, fasta_paths: List[Path], out_prefix: Path, threads: int, stderr_path: Path) -> Path:
    cmd = [mash_cmd, "sketch", "-p", str(threads), "-o", str(out_prefix)] + [str(p) for p in fasta_paths]
    run_checked(cmd, stderr_path=stderr_path)
    msh = out_prefix.with_suffix(".msh")
    if not msh.exists():
        raise RuntimeError(f"mash sketch did not create expected file: {msh}")
    return msh

def parse_mash_dist_output(text: str) -> List[Dict[str, str]]:
    rows = []
    for line in text.strip().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        rows.append({
            "ref": parts[0],
            "query": parts[1],
            "distance": parts[2],
            "p_value": parts[3],
            "shared_hashes": parts[4],
        })
    return rows

def mash_dist(mash_cmd: str, a: Path, b: Path, threads: int) -> List[Dict[str, str]]:
    cmd = [mash_cmd, "dist", "-p", str(threads), str(a), str(b)]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"mash dist failed: {' '.join(cmd)}\n{proc.stderr}")
    return parse_mash_dist_output(proc.stdout)

def build_components(nodes: Iterable[str], edges: Iterable[Tuple[str, str]]) -> List[List[str]]:
    graph: Dict[str, Set[str]] = {n: set() for n in nodes}
    for a, b in edges:
        graph.setdefault(a, set()).add(b)
        graph.setdefault(b, set()).add(a)
    seen = set()
    components = []
    for node in sorted(graph):
        if node in seen:
            continue
        stack = [node]
        comp = []
        seen.add(node)
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for nxt in sorted(graph.get(cur, set())):
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        components.append(sorted(comp))
    return components

def write_svg_network(nodes: Dict[str, int], edges: Dict[Tuple[str, str], int], out_svg: Path) -> None:
   if not nodes:
        out_svg.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="800" height="400"><text x="20" y="40">No replicon network to draw.</text></svg>')
        return
    width, height = 1200, 1200
    cx, cy = width / 2, height / 2
    radius = 420
    ordered = sorted(nodes)
    positions = {}
    for i, node in enumerate(ordered):
        angle = 2 * math.pi * i / max(len(ordered), 1)
        positions[node] = (cx + radius * math.cos(angle), cy + radius * math.sin(angle))
    max_count = max(nodes.values()) if nodes else 1
    max_edge = max(edges.values()) if edges else 1

    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">']
    parts.append('<rect width="100%" height="100%" fill="white"/>')
    parts.append('<text x="30" y="40" font-size="24" font-family="Arial">Replicon-type co-occurrence network from Mash-kept plasmids</text>')
    for (a, b), weight in sorted(edges.items()):
        x1, y1 = positions[a]
        x2, y2 = positions[b]
        stroke_w = 1 + 6 * (weight / max_edge)
        parts.append(
            f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="black" stroke-opacity="0.35" stroke-width="{stroke_w:.2f}"/>'
        )
    for node, count in sorted(nodes.items()):
        x, y = positions[node]
        r = 12 + 28 * (count / max_count)
        parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" fill="lightgray" stroke="black" stroke-width="1.5"/>')
        parts.append(f'<text x="{x:.1f}" y="{y + r + 18:.1f}" text-anchor="middle" font-size="16" font-family="Arial">{node} ({count})</text>')
    parts.append('</svg>')
    out_svg.write_text("\n".join(parts))

def write_dot(nodes: Dict[str, int], edges: Dict[Tuple[str, str], int], out_dot: Path) -> None:
    with open(out_dot, "w") as handle:
        handle.write("graph replicon_network {\n")
        handle.write('  graph [overlap=false, splines=true];\n')
        for node, count in sorted(nodes.items()):
            handle.write(f'  "{node}" [label="{node}\\ncount={count}"];\n')
        for (a, b), weight in sorted(edges.items()):
            handle.write(f'  "{a}" -- "{b}" [label="{weight}"];\n')
        handle.write("}\n")

def exported_value_from_row(row: Dict[str, str]) -> str:
    return (
        row.get("plasmid_fasta_exported")
        or row.get("plasmid_exported")
        or ""
    ).strip().lower()

def exported_path_from_row(row: Dict[str, str]) -> str:
    return (
        row.get("plasmid_fasta_path")
        or row.get("plasmid_export_path")
        or ""
    ).strip()

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Sketch each plasmid group with Mash, keep plasmids having at least one same-group partner "
                    "within the Mash threshold, run pling on kept plasmids, then build a replicon-type network."
    )
    parser.add_argument("--annotation-glob", default="annotation_*")
    parser.add_argument("--output-dir-name", default="mash_then_pling")
    parser.add_argument("--mash-threshold", type=float, default=0.005)
    parser.add_argument("--pling-containment", type=float, default=0.3)
    parser.add_argument("--mash-cmd", default="mash")
    parser.add_argument("--pling-cmd", default="pling")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--link-mode", choices=["symlink", "copy"], default="symlink")
    parser.add_argument("--skip-pling", action="store_true")
    args = parser.parse_args()
   if shutil.which(args.mash_cmd) is None:
        print(f"ERROR: mash command not found: {args.mash_cmd}", file=sys.stderr)
        return 1
    if not args.skip_pling and shutil.which(args.pling_cmd) is None:
        print(f"ERROR: pling command not found: {args.pling_cmd}", file=sys.stderr)
        return 1

    annotation_dirs = [p for p in sorted(Path(".").glob(args.annotation_glob)) if p.is_dir()]
    if not annotation_dirs:
        print(f"ERROR: no annotation directories matched {args.annotation_glob}", file=sys.stderr)
        return 1

    global_group_rows: List[Dict[str, str]] = []
    plasmid_to_replicons_global: Dict[str, Set[str]] = defaultdict(set)
    plasmid_to_annotations: Dict[str, Set[str]] = defaultdict(set)

    for annotation_dir in annotation_dirs:
        try:
            summary_file = choose_summary_file(annotation_dir)
        except FileNotFoundError:
            print(f"[WARN] skipping {annotation_dir}: no summary file", file=sys.stderr)
            continue

        outroot = annotation_dir / args.output_dir_name
        outroot.mkdir(parents=True, exist_ok=True)

        groups: Dict[str, Dict[str, Dict[str, str]]] = defaultdict(dict)
        plasmid_to_replicons_local: Dict[str, Set[str]] = defaultdict(set)

        with open(summary_file, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"sample", "replicon_type"}
            missing = required.difference(reader.fieldnames or [])
            if missing:
                raise RuntimeError(f"{summary_file} missing columns: {', '.join(sorted(missing))}")

            for row in reader:
                if exported_value_from_row(row) != "yes":
                    continue

                fasta_path_str = exported_path_from_row(row)
                if not fasta_path_str:
                    continue
                fasta_path = Path(fasta_path_str)
                if not fasta_path.exists():
                    continue

                replicons = split_replicons((row.get("replicon_type") or "").strip())
                if not replicons:
                    replicons = ["untyped"]

                fasta_key = str(fasta_path.resolve())
                for rep in replicons:
                    plasmid_to_replicons_local[fasta_key].add(rep)
                    plasmid_to_replicons_global[fasta_key].add(rep)
                    plasmid_to_annotations[fasta_key].add(str(annotation_dir))
                    groups[rep][fasta_key] = {
                        "sample": (row.get("sample") or "").strip(),
                        "fasta_path": fasta_key,
                    }

        kept_plasmids_local: Set[str] = set()
        summary_rows: List[Dict[str, str]] = []

        for replicon, plasmid_map in sorted(groups.items()):
            plasmids = sorted([Path(v["fasta_path"]) for v in plasmid_map.values()], key=lambda p: p.name)
            if not plasmids:
                continue

            group_dir = outroot / sanitize(replicon)
            all_dir = group_dir / "all_plasmids"
            ref_dir = group_dir / "reference"
           keep_dir = group_dir / f"mash_within_{args.mash_threshold}"
            fail_dir = group_dir / f"mash_no_partner_{args.mash_threshold}"
            mash_dir = group_dir / "mash"
            pling_dir = group_dir / "pling"

            for d in (all_dir, ref_dir, keep_dir, fail_dir, mash_dir):
                d.mkdir(parents=True, exist_ok=True)

            lengths = {p: read_fasta_length(p) for p in plasmids}
            reference = min(plasmids, key=lambda p: (lengths[p], p.name))

            for p in plasmids:
                link_or_copy(p, all_dir / p.name, args.link_mode)
            link_or_copy(reference, ref_dir / reference.name, args.link_mode)

            group_sketch = mash_sketch(
                args.mash_cmd,
                plasmids,
                mash_dir / "group_sketch",
                args.threads,
                mash_dir / "group_sketch.stderr.log",
            )
            ref_sketch = mash_sketch(
                args.mash_cmd,
                [reference],
                mash_dir / "reference_sketch",
                args.threads,
                mash_dir / "reference_sketch.stderr.log",
            )

            pairwise_rows = mash_dist(args.mash_cmd, group_sketch, group_sketch, args.threads)
            ref_rows = mash_dist(args.mash_cmd, ref_sketch, group_sketch, args.threads)

            keepers: Set[Path] = set()
            pairwise_tsv = mash_dir / "pairwise_mash.tsv"
            with open(pairwise_tsv, "w", newline="") as out:
                out.write("replicon_group\tplasmid_a\tplasmid_b\tlength_a\tlength_b\treference_plasmid\tmash_distance\tshared_hashes\twithin_threshold\n")
                seen_pairs = set()
                for row in pairwise_rows:
                    a = Path(row["ref"])
                    b = Path(row["query"])
                    if a == b:
                        continue
                    pair_key = tuple(sorted([str(a), str(b)]))
                    if pair_key in seen_pairs:
                        continue
                    seen_pairs.add(pair_key)
                    dist = float(row["distance"])
                    within = "yes" if dist <= args.mash_threshold else "no"
                    out.write(
                        f"{replicon}\t{a}\t{b}\t{lengths.get(a, '')}\t{lengths.get(b, '')}\t{reference}\t{dist}\t{row['shared_hashes']}\t{within}\n"
                    )
                    if within == "yes":
                        keepers.add(a)
                        keepers.add(b)

            ref_vs_all_tsv = mash_dir / "reference_vs_all_mash.tsv"
            with open(ref_vs_all_tsv, "w", newline="") as out:
                out.write("replicon_group\treference_plasmid\tquery_plasmid\treference_length\tquery_length\tmash_distance\tshared_hashes\twithin_threshold\n")
                for row in ref_rows:
                    q = Path(row["query"])
                    dist = float(row["distance"])
                    within = "yes" if dist <= args.mash_threshold else "no"
                    out.write(
                        f"{replicon}\t{reference}\t{q}\t{lengths.get(reference, '')}\t{lengths.get(q, '')}\t{dist}\t{row['shared_hashes']}\t{within}\n"
                    )

            for p in plasmids:
dest = (keep_dir if p in keepers else fail_dir) / p.name
                link_or_copy(p, dest, args.link_mode)

            kept_plasmids_local.update(str(p.resolve()) for p in keepers)

            if not args.skip_pling and len(keepers) >= 2:
                input_list = group_dir / "pling_input.txt"
                with open(input_list, "w") as out:
                    for p in sorted(keepers, key=lambda x: x.name):
                        out.write(str((keep_dir / p.name).resolve()) + "\n")
                pling_stdout = group_dir / "pling.stdout.log"
                pling_stderr = group_dir / "pling.stderr.log"
                with open(pling_stdout, "w") as so, open(pling_stderr, "w") as se:
                    proc = subprocess.run(
                        [
                            args.pling_cmd,
                            str(input_list),
                            str(pling_dir),
                            "align",
                            "--cores", str(args.threads),
                            "--containment_distance", str(args.pling_containment),
                        ],
                        stdout=so, stderr=se, check=False, text=True
                    )
                if proc.returncode == 0:
                    pling_status = "ok"
                    pling_message = ""
                else:
                    pling_status = "failed"
                    pling_message = str(pling_stderr)
            elif args.skip_pling:
                pling_status = "skipped_by_user"
                pling_message = ""
            else:
                pling_status = "not_enough_plasmids_after_mash"
                pling_message = ""

            row = {
                "annotation_dir": str(annotation_dir),
                "summary_file": str(summary_file),
                "replicon_group": replicon,
                "n_total_plasmids": str(len(plasmids)),
                "reference_plasmid": str(reference),
                "reference_length_bp": str(lengths[reference]),
                "n_within_mash_threshold": str(len(keepers)),
                "n_without_partner_within_mash_threshold": str(len(plasmids) - len(keepers)),
                "mash_threshold": str(args.mash_threshold),
                "pling_containment_threshold": str(args.pling_containment),
                "pling_status": pling_status,
                "pling_message": pling_message,
                "group_output_dir": str(group_dir),
            }
            summary_rows.append(row)
            global_group_rows.append(row)

        fields = [
            "annotation_dir", "summary_file", "replicon_group", "n_total_plasmids",
           "reference_plasmid", "reference_length_bp", "n_within_mash_threshold",
            "n_without_partner_within_mash_threshold", "mash_threshold",
            "pling_containment_threshold", "pling_status", "pling_message",
            "group_output_dir"
        ]
        with open(outroot / "group_summary.tsv", "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(summary_rows)

        kept_replicon_nodes: Counter = Counter()
        kept_replicon_edges: Counter = Counter()
        for plasmid_path, reps in plasmid_to_replicons_local.items():
            if plasmid_path not in kept_plasmids_local:
                continue
            unique_reps = sorted(set(reps))
            for rep in unique_reps:
                kept_replicon_nodes[rep] += 1
            for a, b in combinations(unique_reps, 2):
                kept_replicon_edges[(a, b)] += 1

        network_dir = outroot / "replicon_network"
        network_dir.mkdir(parents=True, exist_ok=True)
        with open(network_dir / "replicon_nodes.tsv", "w", newline="") as handle:
            handle.write("replicon_type\tkept_plasmid_count\n")
            for rep, count in sorted(kept_replicon_nodes.items()):
                handle.write(f"{rep}\t{count}\n")
        with open(network_dir / "replicon_edges.tsv", "w", newline="") as handle:
            handle.write("replicon_a\treplicon_b\tshared_kept_plasmid_count\n")
            for (a, b), count in sorted(kept_replicon_edges.items()):
                handle.write(f"{a}\t{b}\t{count}\n")

        components = build_components(kept_replicon_nodes.keys(), kept_replicon_edges.keys())
        with open(network_dir / "replicon_clusters.tsv", "w", newline="") as handle:
            handle.write("cluster_id\treplicon_type\n")
            for idx, comp in enumerate(components, start=1):
                for rep in comp:
                    handle.write(f"cluster_{idx}\t{rep}\n")

        write_svg_network(dict(kept_replicon_nodes), dict(kept_replicon_edges), network_dir / "replicon_network.svg")
        write_dot(dict(kept_replicon_nodes), dict(kept_replicon_edges), network_dir / "replicon_network.dot")

    if not global_group_rows:
        print("No groups were processed.", file=sys.stderr)
        return 1

    with open("all_annotations_mash_pling_group_summary.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(global_group_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(global_group_rows)

    global_nodes: Counter = Counter()
    global_edges: Counter = Counter()
    kept_paths = set()
    for row in global_group_rows:
        group_dir = Path(row["group_output_dir"])
        keep_dir = group_dir / f"mash_within_{args.mash_threshold}"
        if keep_dir.exists():
            for fasta in keep_dir.iterdir():
                if fasta.is_file() or fasta.is_symlink():
                    kept_paths.add(str(fasta.resolve()))

    for plasmid_path in kept_paths:
        reps = sorted(set(plasmid_to_replicons_global.get(plasmid_path, [])))
        for rep in reps:
            global_nodes[rep] += 1
        for a, b in combinations(reps, 2):
            global_edges[(a, b)] += 1

    final_dir = Path("all_annotations_replicon_network")
 final_dir.mkdir(parents=True, exist_ok=True)
    with open(final_dir / "replicon_nodes.tsv", "w", newline="") as handle:
        handle.write("replicon_type\tkept_plasmid_count\n")
        for rep, count in sorted(global_nodes.items()):
            handle.write(f"{rep}\t{count}\n")
    with open(final_dir / "replicon_edges.tsv", "w", newline="") as handle:
        handle.write("replicon_a\treplicon_b\tshared_kept_plasmid_count\n")
        for (a, b), count in sorted(global_edges.items()):
            handle.write(f"{a}\t{b}\t{count}\n")

    components = build_components(global_nodes.keys(), global_edges.keys())
    with open(final_dir / "replicon_clusters.tsv", "w", newline="") as handle:
        handle.write("cluster_id\treplicon_type\n")
        for idx, comp in enumerate(components, start=1):
            for rep in comp:
                handle.write(f"cluster_{idx}\t{rep}\n")

    write_svg_network(dict(global_nodes), dict(global_edges), final_dir / "replicon_network.svg")
    write_dot(dict(global_nodes), dict(global_edges), final_dir / "replicon_network.dot")

    print("Wrote all_annotations_mash_pling_group_summary.tsv")
    print(f"Wrote final replicon network to {final_dir}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())


