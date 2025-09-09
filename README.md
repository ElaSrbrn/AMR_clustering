# CRE pipeline
# plasmid-outbreak-pipeline

Bioinformatic pipeline for Nanopore plasmid and chromosomal clustering analysis: basecalling → demux → filter_qc → assembly → polishing → plasmid & AMR annotation → plasmid clustering and MLST typing (on downsampled coverage files).

Make sure to update the Library prep kit, depending on the run. For MinION run, some concatenated files are not consecutive, so deduplication and concatenation might have to be done manually.

```

