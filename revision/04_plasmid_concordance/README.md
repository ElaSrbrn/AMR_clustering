# Plasmid contig stability pipeline

Compares the carbapenemase-carrying plasmid contigs of Flye and Hybracter assemblies (10/20/40/60x)
with the Autocycler 60x reference (with Plassembler), classifies how each contig deviates, selects
plasmids for in-depth analysis and synteny figures, and prepares the input for the Mash/DCJ clustering.

Plasmid contig stability = structural concordance of the carbapenemase-carrying contig with the
reference plasmid: |length difference| <= 0.1%, >= 99% of the reference aligned in one collinear
block (no breakpoint), identity >= 99.9%. Thresholds are in `config.sh`.

## Location

Put this folder at `/ictstr01/groups/hpc/urban_lab/projects/plasmid_project/plasmid_contig_stability`
and run every command from inside it.

```
plasmid_contig_stability/
  config.sh                   paths, environments, databases, thresholds  (fill PLACEHOLDER lines)
  ids_s1a_60x.txt             98 carbapenemase-positive isolates with a 60x reference
  overrides.tsv               optional explicit assembly paths (id, assembler, depth, path)
  manual_include_synteny.tsv  plasmids to add to the synteny proposal (e.g. Mash/DCJ cluster medoids)
  lib/common.py               shared Python helpers (standard library only)
  scripts/                    numbered steps
  work/                       intermediate files (manifest, annotation, alignments, in-depth folders)
  results/                    tables for the paper and the response letter
  selection/                  YOUR edited copies of the proposals (never overwritten by the scripts)
  synteny_packages/           one folder per plasmid for the collaborator
  clustering_input/           FASTA sets for Mash/DCJ
  logs/                       SLURM logs
```

## Before the first run

1. Wait until all Autocycler v4 and Hybracter reruns report `OK`, so the references and test assemblies are final.
2. In `config.sh`, fill in:
   - `FLYE_GLOB`: where the Medaka-polished Flye assemblies are. `bash scripts/00_find_flye.sh 6817` lists candidate files for one isolate.
   - `BAKTA_DB`, `ISFINDER_DB`: as used in the EMBL project.
   - `IS_MIN_ID`, `IS_MIN_COV`: replace with the EMBL ISfinder thresholds.
3. For isolates with more than one read set (6901: use FBC05034 at every depth), add explicit paths to `overrides.tsv` if the manifest reports `MULTIPLE`.

## Run order

```bash
source config.sh

# 00  manifest of all assemblies; fix every MISSING / MULTIPLE / EMPTY row before continuing
bash scripts/00_make_manifest.sh

# 01  AMRFinderPlus + MOB-typer on every assembly, one task per isolate (resubmittable, finished ones are skipped)
mkdir -p logs
sbatch --array=1 scripts/01_annotate.sbatch
sbatch --array=2-$(wc -l < ids_s1a_60x.txt)%15 scripts/01_annotate.sbatch
grep -h "^STATUS" logs/01_annotate_*.out          # every isolate should show n/n

# 02-04  reference plasmids, pairing, alignment + classification, summaries
sbatch scripts/02_03_collect_compare.sbatch

# 05  proposals for in-depth analysis and synteny; review, copy, edit 'include'
$PY scripts/05_select.py
mkdir -p selection
cp results/proposal_indepth.tsv selection/indepth.tsv
cp results/proposal_synteny.tsv selection/synteny.tsv
#     edit selection/*.tsv: include = yes / no

# 06  Bakta + ISfinder + self-repeats for the selected plasmids
$PY scripts/06_prepare_indepth.py
sbatch --array=1-$(wc -l < work/indepth_list.txt)%10 scripts/06_indepth.sbatch

# 07  IS elements and repeats vs breakpoints
$PY scripts/07_breakpoints_vs_is.py

# 08  synteny packages for the collaborator
$PY scripts/08_package_synteny.py
tar -czf synteny_packages.tar.gz synteny_packages

# 09  FASTA sets for the Mash/DCJ (Pling) clustering: 60x reference vs Flye 40x (add more sets as arguments)
$PY scripts/09_export_for_clustering.py flye_40x
```

## Main outputs

| File | Content |
|---|---|
| `results/reference_carbapenemase_loci.tsv` | every carbapenemase gene in the references, plasmid or chromosome (basis for S1a) |
| `results/reference_plasmids.tsv` | reference plasmids entering the analysis: genes, length, circularity, replicon, mobility |
| `results/stability_long.tsv` | one row per plasmid x assembler x depth: class and all metrics |
| `results/stability_rates.tsv` | proportion stable per assembler and depth (strict = circular; structural = circularity not required) |
| `results/summary_by_depth.tsv` | counts of each class per assembler and depth |
| `results/stability_matrix.tsv` | plasmid x assembler_depth -> class (heatmap input) |
| `results/breakpoints.tsv` | breakpoints in reference coordinates; `size_bp` > 0 = reference sequence missing in the test contig, < 0 = extra sequence in the test contig |
| `results/breakpoints_annotated.tsv`, `results/breakpoint_is_summary.tsv` | distance of breakpoints to IS elements/repeats, observed vs expected proportion |
| `results/is_elements_reference.tsv`, `results/is_elements_test_contigs.tsv` | ISfinder hits |

## Classes (first matching rule wins)

`gene_absent`, `on_chromosome_contig`, `concordant`, `concordant_not_circular`, `rearranged`,
`expanded_duplication` (duplicated reference sequence, e.g. overlapping contig ends), `chimeric`
(> 5% of the contig does not align), `internal_deletion`, `internal_insertion`, `fragmented`
(the rest of the plasmid is on other contigs), `truncated`, `length_mismatch`, `low_identity`.
Full definitions at the top of `scripts/03_compare.py`.

## Notes and limitations

- The test contig compared is the one carrying the carbapenemase gene. If several contigs carry it, the one covering most of the reference is used; `n_gene_contigs` records how many there were.
- A correct circular contig with a different start position, or on the opposite strand, is not counted as a breakpoint.
- Plasmids sharing backbone with another plasmid in the same isolate can inflate `n_fragment_contigs`; check fragmented calls in multi-plasmid isolates by hand.
- Breakpoint/IS association: the two ends of one breakpoint are not independent. Compare observed with `expected_prop_random`, and test per breakpoint or per plasmid.
- Flye is not fully deterministic when multithreaded, so borderline calls at low depth can differ between runs. Report proportions across plasmids.
- Tested on synthetic assemblies covering every class (truncation, fragmentation, rotation on both strands, overlapping ends, insertion, deletion, inversion, gene loss). Not yet run on the real data.
