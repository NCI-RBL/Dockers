### Prepare directories

Create a working directory and subdirectories:

```bash
mkdir DTEG
cd DTEG
mkdir FASTQ
mkdir HTSeq
```


Create a subdirectory for each sample. Example:

```bash
mkdir FASTQ/HEK_WT1
mkdir FASTQ/HEK_WT2
mkdir FASTQ/HEK_ThumD_KO1
mkdir FASTQ/HEK_ThumD_KO2
mkdir FASTQ/HEK293ThumpD1KO1RNA
mkdir FASTQ/HEK293ThumpD1KO2RNA
mkdir FASTQ/HEK293WT1RNA
mkdir FASTQ/HEK293WT2RNA
```


### Index

You can link the pre-made index folder for hg19 in the working directory:

```bash
ln -s /mnt/rnabl-work/Guiblet/CCBRRBL8/NextFlow ./
```

You can create a new index as follow:
TBA


### Scripts

Download the following files in the working directory:

- [DTEG_RBL.nf](https://github.com/RBL-NCI/Dockers/blob/main/workflows/DTEG/DTEG_RBL.nf)
- [DTEG_RBL.sh](https://github.com/RBL-NCI/Dockers/blob/main/workflows/DTEG/DTEG_RBL.sh)
- [nextflow.config](https://github.com/RBL-NCI/Dockers/blob/main/workflows/DTEG/nextflow.config)
