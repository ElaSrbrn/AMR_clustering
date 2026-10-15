# CRE pipeline

Bioinformatic pipeline for nanopore sequencing of carbapenem-resistant *Enterobacterales*, with a focus on plasmid-encoded carbapenem resistance and on how much sequencing depth each analysis needs.

This repository holds the code for the manuscript "Verification of nanopore sequencing technology for clinical carbapenem-resistant Enterobacterales surveillance" (preprint: bioRxiv, doi.org/10.64898/2026.08.04.742736; revised version under review at the Journal of Medical Microbiology). Filtered sequencing reads are deposited at ENA under project `PRJEB123044`; sample and run accessions are listed in Supplementary Table S1.

The repository has two parts:

- **Top-level scripts (00–09):** the core pipeline, from basecalling to plasmid clustering, as used in the original submission and unchanged in the revision.
- **`revision/`:** the analyses added during peer review (reference assemblies, Hybracter, structural plasmid concordance, read-based typing, cgMLST, clustering at 40×, statistics), with their result tables. See `revision/README.md`.

Workflow: basecalling → read filtering → nested downsampling → assembly → polishing → AMR and plasmid annotation → MLST → plasmid clustering, with typing and plasmid analyses repeated at each depth. All depth analyses use the 100 carbapenemase-producing isolates (Supplementary Table S1a); the 10 carbapenem-non-susceptible isolates (S1b) were characterised from their 40× read sets only.

## Core pipeline: running order

| Script | Does | Needs |
|---|---|---|
| `00_basecalling.sh` | Dorado basecalling, POD5 to BAM | POD5 from the sequencer |
| `01_convert_bam_tofastq.sh` | BAM to FASTQ | output of 00 |
| `02_chopper.sh` | filters to ≥1000 bp and mean Q≥9, drops duplicate read IDs | `./fastq/*barcode*.fastq.gz` |
| `03_downsampling_rasusa_nested.sh` | nested downsampling to each target depth | output of 02 |
| `04_assembly.sh <depth>` | Flye assembly, one depth per call | `./<depth>x/` |
| `05_polishing.sh <depth>` | Medaka polishing | output of 04 |
| `06_amr_plasmid_annotation.sh` | AMRFinderPlus, Bakta, MOB-suite | polished assemblies |
| `07_mlst.sh` | MLST typing at each depth | polished assemblies |
| `08_clustering.sh` | Mash filtering then pling clustering | `annotation_*` from 06 |
| `09_bart_amr_mlst.sh` | stand-alone wrapper for read-based MLST and AMR typing with bart | read folders per depth |
| `build_carbapenemase_summary.py` | collates carbapenemase calls across depths | output of 06 |

Scripts 04 and 05 take the depth directory as their first argument, for example `sbatch 04_assembly.sh 60x`.

## Environments

```bash
micromamba create -f assembly.yml   # chopper, Flye, Medaka, minimap2, samtools, SeqKit
micromamba create -f amr.yml        # AMRFinderPlus, Bakta
micromamba create -f mlst.yml       # mlst, BLAST
micromamba create -f mobsuite.yml   # MOB-suite, Mash
micromamba create -f pling.yml      # pling
```

The revision analyses additionally need Autocycler, Hybracter, chewBBACA and bart; versions are listed below.

## Versions used (as in the revised Methods)

Dorado v1.1.1 (model dna_r10.4.1_e8.2_400bps_sup@v5.0.0), Chopper v0.12, SeqKit v2.10.0, Rasusa v3.0.0 (seed 42), Flye v2.9.6, Medaka v2.1.1, Minimap2 v2.31, CheckM2 v1.1.0, Pathogenwatch v2.3.1, mlst v2.23.0, AMRFinderPlus v4.2.7 (database 2026-05-15.1), MOB-suite v3.1.8, Mash v2.3, pling v2.0.1. Revision analyses: Autocycler v0.8.0, Plassembler v1.8.5, Hybracter v0.14.1, chewBBACA v3.5.4, bart v0.1.2 with KMA v1.6.17, Bakta v1.11.4, Python v3.13 with statsmodels v0.15.0 and SciPy v1.18.1.

## Key parameters

Reads were filtered to ≥1000 bp and mean quality ≥Q9. Downsampling used a genome size of 5 Mb to target nominal depths of 100×, 80×, 60×, 40×, 20× and 10×, nested so that each lower depth is a subset of the one above. Actual depths relative to each isolate's reference assembly are given in Supplementary Table S6. Assemblies used Flye `--nano-hq` and Medaka in bacterial mode.

For plasmid clustering, plasmids of the same replicon type were retained if within a Mash distance of 0.005 of at least one other plasmid, then clustered with pling at a containment distance of 0.3 and a DCJ-Indel threshold of 4.

## Notes from doing this

We duplexed two barcodes per sample, which is not strictly necessary. If you do the same, concatenate the filtered reads and add a deduplication step before assembly. The library prep kit needs updating in the basecalling script depending on the run.

For downsampling we used a conservative genome size of 5 Mb because it is easy to apply across species. Actual depth was 0.79 to 1.26 times the nominal depth, lowest in *Klebsiella* species with genomes above 5 Mb (Supplementary Table S6).

The pipeline is written for a SLURM cluster. Partition and queue names in the SBATCH headers, and the absolute paths in `revision/`, are specific to our system and need changing.

## Data

Supplementary Tables 4 and 5 summarise the Mash distance and DCJ-Indel output produced by `08_clustering.sh`, and Figure 4 was drawn from the Flye 60× output. The statistics in the revised manuscript can be reproduced from the Supplementary Tables workbook with `revision/10_statistics/manuscript_statistics.py`.

## Citation

If you use this code, please cite the manuscript (citation to be added on publication) or the preprint: doi.org/10.64898/2026.08.04.742736.
