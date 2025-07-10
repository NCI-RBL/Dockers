### Input files

You will need:

- A STAR index directory (example: STAR_hg38)
- A genome annotation (example: gencode.v44.annotation.gtf)
- FASTQ files (we will assume paired files for this tutorial)


### Index

If you do not have a STAR index you can create one:

```bash
STAR \
  --runThreadN 64 \
  --runMode genomeGenerate \
  --genomeDir STAR_hg38/ \
  --genomeFastaFiles references/hg38.fa \
  --sjdbGTFfile references/gencode.v44.annotation.gtf
```

