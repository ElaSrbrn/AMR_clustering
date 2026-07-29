# CRE pipeline

Bioinformatic pipeline for nanopore sequencing of carbapenem-resistant *Enterobacterales*, with a focus on plasmid-encoded carbapenem resistance and on how much sequencing depth each analysis actually needs.
This pipeline was used for the manuscript on "Verification of nanopore sequencing technology for clinical carbapenem-resistant Enterobacterales surveillance" submitted to BioRxiv and available here: XXXX

Workflow: basecalling → read filtering → nested downsampling → assembly → polishing → AMR and plasmid annotation → MLST → plasmid clustering, with the typing and plasmid analyses repeated at each depth.

This repository holds the code for the manuscript below. Filtered sequencing reads are deposited at ENA under accession `PRJEBXXXXXX`.

## Running order

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
| `build_carbapenemase_summary.py` | collates carbapenemase calls across depths | output of 06 |
| `st_depth_analysis.Rmd` | ST concordance versus depth | per-depth MLST tables from 07 |
| `plasmid_stability_depth_analysis.Rmd` | plasmid reconstruction stability versus depth | plasmid size table from 06 |

Scripts 04 and 05 take the coverage directory as their first argument, for example `sbatch 04_assembly.sh 60x`.

## Environments

Four micromamba environments, one YAML each:

```bash
micromamba create -f assembly.yml   # chopper, Flye, Medaka, minimap2, samtools, SeqKit
micromamba create -f amr.yml        # AMRFinderPlus, Bakta
micromamba create -f mlst.yml       # mlst, BLAST
micromamba create -f mobsuite.yml   # MOB-suite, Mash
micromamba create -f pling.yml      # pling
```

The R analyses need R v4.3.1 with `dplyr` v1.1.4, `tidyr`, `geepack` v1.3.12, `logistf` v1.26.1, `sandwich`, `lmtest`, `ggplot2` v3.5.2 and `pheatmap` v1.0.13.

## Versions used

Dorado v1.1.1, chopper v0.12, Rasusa v3.0.0 (seed 42), Flye v2.9.6, Medaka v2.1.1, minimap2, CheckM2 v1.1.0, SeqKit v2.10.1, mlst v2.32.0, AMRFinderPlus v4.0.23 (database v4.0.1), MOB-suite v3.1.8, Mash v2.3, pling v1.0.1, R v4.3.1.

## Key parameters

Reads were filtered to ≥1000 bp and mean quality ≥Q9. Downsampling used a genome size estimate of 5 Mb to target depths of 100×, 80×, 60×, 40×, 20× and 10×, nested so that each lower depth is a subset of the one above. Assemblies used Flye `--nano-hq` and Medaka in bacterial mode.

For plasmid clustering, plasmids of the same replicon type were retained if within a Mash distance of 0.005 of at least one other plasmid, then clustered with pling at a containment distance of 0.3 and a DCJ-Indel threshold of 4.

## Notes from doing this

We duplexed two barcodes per sample, which is not strictly necessary. If you do the same, concatenate the filtered reads and add a deduplication step before assembly. The library prep kit needs updating in the basecalling script depending on the run.

For downsampling we used a conservative genome size estimate of 5 Mb because it is easy to apply across species. In practice the estimated sequencing depth did not always match the median chromosomal coverage, and coverage generally ran lower than the nominal depth, partly because of high-copy-number plasmids. Using 5.2 Mb would make 10× correspond more closely to 10× chromosomal coverage.

The pipeline is written for a SLURM cluster. Partition and queue names in the SBATCH headers are specific to our system and will need changing.

## Data

Supplementary Tables S1, S3 and S6 from the manuscript can be used as input for the R analyses. Otherwise, the necessary data for the MLST analysis is included in the repository as well. 
Supplementary Tables 4 and 5 summarise the Mash distance and DCJ-Indel output produced by 08_clustering.sh, and Figure 4 was drawn from the same output. No separate script is included for those steps.

## Citation

If you use this code, please cite:

> XXXXX

