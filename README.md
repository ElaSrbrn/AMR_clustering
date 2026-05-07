# CRE pipeline
# plasmid-annotation-pipeline

Bioinformatic pipeline for Nanopore plasmid and chromosomal clustering analysis: basecalling  → chopper → downsampling -> assembly → polishing → plasmid & AMR annotation → MLST -> plasmid clustering and MLST typing for different depths.

Make sure to update the Library prep kit, depending on the run. For downsampling, we have decided on a conservative and easy to implement genome size estimate of 5Mb. However, we have seen that sequencing depth calculated on the estimated sequencing depth did not always match the median chromosomal coverage with the coverage mainly being lower than the sequencing depth, also due to high copy number plasmids. An alternative is to use e.g. 5.2Mb as estimate, this would ensure that 10x sequencing depth corresponds more closely to 10x coverage. 



