### RiboFootPrint

Steps of the pipeline:

- Trimm Ribo-seq reads.
- Filter out reads mapping to non-coding RNA and remap (STAR) to transcriptome.
- Aggregates relative postions (range: 0-1) of read starts, split in 5'UTR, CDS, and 3'UTR.

For more questions about this pipeline, you may contact [Colin Wu](mailto:colin.wu2@nih.gov) or [Wilfried Guiblet](mailto:guibletwm@nih.gov).
