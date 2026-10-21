# Analyses added for the revision

This folder holds the code and result tables of the analyses added during peer review of the manuscript (J Med Microbiol). The core pipeline (basecalling, filtering, downsampling, Flye/Medaka assembly, annotation, MLST, Mash/pling at 60×) is in the top-level scripts and was not changed.

All depth analyses cover the 100 carbapenemase-producing isolates of Supplementary Table S1a. Plasmid analyses cover the 98 of them with a 60× read set (6804 and 7418 have none) and their 93 carbapenemase-encoding reference plasmids in 89 isolates. Isolate 7418 reached a maximum of 30×.

`COPY_LOG.tsv` gives the cluster path each file was copied from; `MISSING.txt` lists anything that was not found when the folder was built. Absolute paths in the scripts refer to our cluster:

```
BASE   = /ictstr01/groups/hpc/urban_lab/projects/plasmid_project
NEWCRE = $BASE/stored_resistant_enterobacterales/new_cre
PCS    = $BASE/plasmid_contig_stability
```

Reads, assemblies, alignments and annotation output are not included. Raw data: ENA PRJEB123044.

| Folder | Analysis | Manuscript |
|---|---|---|
| 01_read_sets_by_depth | read sets per isolate and depth; Hybracter sample sheet | Methods |
| 02_autocycler_references | Autocycler 60× reference assemblies with Plassembler | Methods; S3 |
| 03_hybracter | Hybracter at 10×, 20×, 40× and 60× | Methods; S3; Table 1 |
| 04_plasmid_concordance | structural concordance, actual depth, IS annotation, read support, tandem copy number, context synteny, PLSDB comparison, Hybracter gene loss, clustering of the Autocycler references | Results; Table 1; S3, S6, S9, S10b, S11, S12 |
| 05_synteny | input description for the synteny figures | Figure 5 |
| 06_read_typing_bart | carbapenemase and MLST typing from reads | Results; Table 1; S8a–c |
| 07_mlst_by_depth | assembly-based MLST at every depth | Results; Table 1; S8d |
| 08_cgmlst | cgMLST reproducibility and clusters by depth | Results; Table 1; S7 |
| 09_plasmid_clustering_40x | Mash and pling on the Flye 40× plasmids, compared with 60× | Results; S10a |
| 10_statistics | all statistics reported in the manuscript | Results; Table 1 |
| 11_figures | figure code | Figures |

## 01_read_sets_by_depth

`build_depth_inputs.sh` builds symlinked read folders per depth for the S1 isolates and writes the Hybracter sample sheet `hybracter_depth_series.csv` and `s1_min_chrom.tsv`, the per-species minimum chromosome length passed to Hybracter with `-c`. Reads were downsampled by `03_downsampling_rasusa_nested.sh` (top level).

## 02_autocycler_references

Autocycler v0.8.0 consensus assemblies from the 60× reads. For each isolate, two read subsets (`autocycler subsample --count 2`) were each assembled with Flye v2.9.6, Raven v1.8.3, miniasm v0.3/Minipolish v0.2.1, metaMDBG v1.4 and myloasm v0.7.0 (`carblot_pos/`, `carblot_neg/`). `rerun_autocycler_plassembler_v4.sh` regenerates the same subsamples, adds Plassembler v1.8.5 assemblies (PLSDB 2023_11_03_v2) and rebuilds the consensus (`compress`, `cluster`, `trim`, `resolve`, `combine`). `reference_check.sh` compares the references with and without Plassembler.

```bash
MODE=check bash rerun_autocycler_plassembler_v4.sh
sbatch --array=1-$(wc -l < isolate_dirs.txt)%20 rerun_autocycler_plassembler_v4.sh
```

## 03_hybracter

Hybracter v0.14.1 (`hybracter long-single`) for every carbapenemase-producing isolate and depth, with `--skip_qc` (reads already filtered), the species-specific `-c` from `s1_min_chrom.tsv` and default Medaka polishing. When Hybracter completes the chromosome, the final assembly contains the chromosome and the Plassembler plasmids; otherwise it contains all polished contigs.

```bash
sbatch --array=1-$(wc -l < hybracter_depth_series.csv)%20 run_hybracter_medaka_array.sh
```

## 04_plasmid_concordance

Numbered pipeline steps in `scripts/`, paths and thresholds in `config.sh`, details in this folder's `README.md`. Results are in `results/`.

| Step | What it does | Main output | Manuscript |
|---|---|---|---|
| 00–01 | manifest of Flye and Hybracter assemblies; AMRFinderPlus and MOB-typer on every assembly | work/ | — |
| 02–04 | reference plasmids; alignment of each carbapenemase-carrying test contig to its reference (Minimap2 asm5); concordance classes; summaries | reference_plasmids.tsv, stability_long.tsv, stability_matrix.tsv, summary_by_depth.tsv | S3a, S3b, Table 1 |
| 06–07 | Bakta v1.11.4 and ISfinder (BLASTn; ≥95% identity, ≥80% coverage) on 59 plasmids; distance of each boundary of a structural difference to IS elements and repeats | is_elements_*.tsv, breakpoints_annotated.tsv, breakpoint_is_summary.tsv | S11a–d |
| 08, 08b | input packages for the synteny figures | ../05_synteny | Figure 5 |
| 10 | actual depth of each read subset (bases / reference length) | actual_depth*.tsv | S6 |
| 12, 12b, 13 | read support and boundary-spanning reads for non-concordant plasmids | read_support.tsv, junction_support.tsv | S9a, S9b |
| 14, 14b–d | nearest PLSDB plasmids (Mash), pairwise alignment to the nearest hit, BioSample metadata | plsdb_*.tsv | S12 |
| 15 | copy number of the 9.3 kb tandem unit of plasmid 7365_2 in 60× reads spanning the array | array_copynumber_7365_2.tsv | S9c |
| 16 | synteny of the carbapenemase gene ±10 kb with the reference in each Flye and Hybracter assembly | context_synteny_*.tsv | S3d |
| 17 | Mash and pling on the Autocycler reference plasmids, compared with the Flye 60× clustering | clustering_autocycler60x/ | S10b |
| 18 | Hybracter assemblies lacking the carbapenemase gene: reference plasmid and gene window in the intermediate Flye assembly, the Plassembler output and the final assembly | hybracter_gene_loss_check.tsv | S3c |

File and column names still say "breakpoint"; the manuscript calls these structural differences and their two ends boundaries. The statistical test of the association between structural differences and IS elements (former step 07b) was not retained in the manuscript and is not included.

## 05_synteny

`SYNTENY_INSTRUCTIONS.txt` describes the input packages built by steps 08 and 08b (reference plasmid and test contigs per depth and assembler, oriented to the reference, with Bakta annotation, IS elements and structural differences). The figures were drawn with pyGenomeViz [version; figure code to be added].

## 06_read_typing_bart

Carbapenemase genes and MLST alleles typed from the 10× and 20× reads with bart v0.1.2 (KMA v1.6.17).

- `build_bart_amr_index.py` builds the AMR index from the AMRFinderPlus database 2026-05-15.1; carbapenemase alleles in families with non-carbapenemase members (such as blaOXA) form a separate family (`check_carbapenemase_families.sh`).
- `run_bart_pipeline.sh` runs `bart_array.sh` and `bart_job.sh`; MLST schemes are mapped in `scheme_map.tsv`.
- `bart_combine.py` merges the outputs, `classify_isolates.py` assigns S1a/S1b, and `compare_reads_vs_assemblies.py` compares the reads with AMRFinderPlus on the highest-depth assembly of the same isolate.
- Calling rule (relative coverage cutoff): per gene family, the top-ranked allele is called if its read coverage is at least 0.2 times the median read coverage of the isolate's core AMR genes. It uses no reference information.
- `results/` holds the deduplicated results (isolate 6901 counted once).

## 07_mlst_by_depth

`all_coverages_mlst_combined.tsv`: mlst v2.23.0 calls on the Medaka-polished Flye assemblies at every depth (from `07_mlst.sh`). Supplementary Table S8d compares each call with the ST in Supplementary Table S1 (40× assembly; 30× for 7418).

## 08_cgmlst

chewBBACA v3.5.4 with cgMLST.org schemes (`schemes.tsv`) on the Medaka-polished Flye assemblies of the 90 carbapenemase-producing isolates with a scheme, at 10×, 20×, 40×, 60×, 80× and 100×. The reference for each isolate is its highest-depth assembly of at least 60× (88 isolates).

```bash
python3 cgmlst_manifest.py      # manifest of assemblies per scheme and depth
bash cgmlst_run.sh prep         # one SLURM job per scheme
bash cgmlst_run.sh call
python3 cgmlst_depth.py         # profiles vs reference, pairwise distances, clusters (--exclude for isolates to leave out)
python3 cgmlst_add_st.py        # adds 7-gene STs
```

The exact arguments are in each script's header. `results/cgmlst_cpe_*` are the final outputs used in the manuscript (carbapenemase-producing isolates only). Linkage thresholds (cgMLST.org complex type distance): *K. pneumoniae* species complex ≤15, *C. freundii* ≤10, *E. coli* ≤10, *K. oxytoca* complex ≤6; none for *E. hormaechei*, *P. mirabilis* and *S. marcescens*.

## 09_plasmid_clustering_40x

`run_annotation_carbapenemase_fixed.slurm` annotates the Flye 40× assemblies; `run_mash_then_pling_by_group_v2.slurm` runs `group_plasmids_mash_then_pling_v2_compat.py` (Mash ≤0.005, then pling v2.0.1 with containment 0.3 and DCJ-Indel ≤4, per replicon group). `compare_pling_depths.py` compares the 40× and 60× assignments pair by pair (`pling_depth_comparison_annotation_60x_vs_annotation_40x.tsv`, S10a).

## 10_statistics

`manuscript_statistics.py` reproduces the statistics in the manuscript from the Supplementary Tables workbook:

- structural concordance per assembler and depth with exact Clopper–Pearson 95% CIs (S3b);
- exact McNemar tests between consecutive depths and between assemblers;
- logistic GEE on log2(depth), exchangeable working correlation clustered by plasmid, robust standard errors, per assembler and jointly with assembler as covariate;
- exact McNemar test of read-based vs assembly-based carbapenemase typing (S8a);
- actual vs nominal depth (S6a).

```bash
python3 manuscript_statistics.py Supplementary_Tables_S1-S12.xlsx statistics_output/
```

Tested with Python 3.13, pandas 2, statsmodels 0.15.0, SciPy 1.18.1 and openpyxl 3.1.


