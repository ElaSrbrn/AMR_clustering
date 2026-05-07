# CRE pipeline

Bioinformatic pipeline for Nanopore plasmid and chromosomal clustering analysis with a focus on plasmid-encoded carbapenem resistance: basecalling  → chopper → downsampling -> assembly → polishing → plasmid & AMR annotation → MLST -> plasmid clustering and MLST typing for different depths.

In our case, we duplexed two barcodes per sample, which is not strictly necessary. If you decide to do so as well, I would recommend concatenating the filtered reads, then adding a deduplication step before assembling them. Make sure to update the Library prep kit, depending on the run for the basecalling algorithm. For downsampling, we have decided on a conservative, easy-to-implement genome size estimate of 5 Mb. However, we have seen that the estimated sequencing depth did not always match the median chromosomal coverage, with coverage generally lower than the estimated sequencing depth, also due to high-copy-number plasmids. An alternative is to use, e.g., 5.2 Mb as an estimate; this would ensure that 10x sequencing depth more closely corresponds to 10x coverage. 



