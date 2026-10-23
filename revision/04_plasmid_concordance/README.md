# Plasmid concordance pipeline

Compares the carbapenemase-carrying plasmid contigs of the Flye and Hybracter assemblies at 10×, 20×, 40× and 60× with the Autocycler 60× reference plasmids, and runs the follow-up analyses on these plasmids (IS elements, read support, context synteny, PLSDB, clustering of the references).

A test contig is structurally concordant with its reference plasmid when the length differs by at most 0.1%, one collinear alignment covers at least 99% of the reference, identity is at least 99.9% and the contig is circular (concordant_not_circular otherwise). All thresholds, paths and environments are set in `config.sh`. Run every step from this folder after `source config.sh`.

```
config.sh                    paths, environments, databases, thresholds
ids_s1a_60x.txt              the 98 carbapenemase-producing isolates with a 60× reference
overrides.tsv                explicit assembly paths where the glob in config.sh is ambiguous (isolate 6901)
manual_include_synteny.tsv   plasmids added by hand to the synteny selection
synteny_groups.tsv           groups of reference plasmids for the comparison figures
lib/common.py                shared helpers (standard library only)
scripts/                     numbered steps
results/                     result tables
```

## Steps

| Step | Script | Output (results/) |
|---|---|---|
| 00 | `00_make_manifest.sh`: one row per isolate, assembler and depth | work/manifest.tsv |
| 01 | `01_annotate.sbatch`: AMRFinderPlus and MOB-typer on every assembly | work/ |
| 02–04 | `02_03_collect_compare.sbatch`: reference plasmids, alignment of each test assembly (Minimap2 asm5), concordance classes, summaries | reference_plasmids.tsv, stability_long.tsv, stability_matrix.tsv, stability_rates.tsv, summary_by_depth.tsv, breakpoints.tsv |
| 05 | `05_select.py`: proposals for in-depth annotation and synteny figures, edited copies in selection/ | proposal_*.tsv |
| 06 | `06_prepare_indepth.py`, `06_indepth.sbatch`: Bakta, ISfinder (BLASTn) and self-repeats per selected plasmid | work/indepth/ |
| 07 | `07_breakpoints_vs_is.py`: distance of each boundary of a structural difference to IS elements and repeats | breakpoints_annotated.tsv, breakpoint_is_summary.tsv, is_elements_*.tsv |
| 08, 08b | `08_package_synteny.py`, `08b_synteny_groups.py`: input folders for the synteny figures | ../synteny_packages/, ../synteny_groups/ |
| 09, 11 | `09_export_for_clustering.py`, `11_compare_pling_depths.py`: plasmid FASTA sets and comparison of clusterings | clustering_input/ |
| 10 | `10_actual_depth.sbatch`: total bases of each read subset divided by the reference length | actual_depth.tsv, actual_depth_summary.tsv |
| 12, 13 | `12_read_support.sbatch`, `13_junction_support.sbatch`: 60× reads mapped to reference and Flye contig; reads spanning each boundary | read_support.tsv, junction_support.tsv |
| 14 | `14_plsdb.sbatch`, `14b_plsdb_align.sbatch`, `14c_plsdb_meta.py`, `14d_plsdb_provenance.py`: nearest PLSDB plasmid (Mash), alignment to it, BioSample metadata | plsdb_*.tsv |
| 15 | `15_array_copynumber.sbatch`: copy number of the tandem unit of plasmid 7365_2 in reads spanning the array | array_copynumber_7365_2*.tsv |
| 16 | `16_context_synteny.sbatch`: carbapenemase gene ±10 kb aligned to every test assembly | context_synteny_*.tsv |
| 17 | `17_autocycler_pling.sbatch`: Mash and pling on the Autocycler reference plasmids, compared with the Flye 60× clustering | ../clustering_autocycler60x/ |
| 18 | `18_hybracter_gene_loss_check.py`: Hybracter assemblies without the carbapenemase gene, checked in the intermediate Flye assembly and the Plassembler output | hybracter_gene_loss_check.tsv |

Concordance classes, in the order the rules are applied: gene_absent, on_chromosome_contig, concordant, concordant_not_circular, rearranged, expanded_duplication, chimeric, internal_deletion, internal_insertion, fragmented, truncated, length_mismatch, low_identity. Definitions are at the top of `scripts/03_compare.py` and in Supplementary Table S3.

File and column names say "breakpoint"; the manuscript calls these structural differences, and their two ends boundaries.
