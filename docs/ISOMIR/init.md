### Prepare directories

Create a working directory and subdirectories:

```bash
mkdir ISOMIR
cd ISOMIR
mkdir fastq_files
```


### Index

You will need an Index directory for mapping. You can download this [Index](https://hpc.nih.gov/~RBL_NCI/Dockers/isomir/index.tar.gz) for hg38.
Extract the downloaded directory :

```bash
tar zxvf index.tar.gz
```


### Scripts

Download the following files in the working directory:

- [isomiR.nf](https://github.com/NCI-RBL/Dockers/blob/main/workflows/isomiR/isomiR.nf)
- [isomiR.sh](https://github.com/NCI-RBL/Dockers/blob/main/workflows/isomiR/isomiR.sh)
- [nextflow.config](https://github.com/NCI-RBL/Dockers/blob/main/workflows/isomiR/nextflow.config)
- [motif-consensus.fa](https://github.com/NCI-RBL/Dockers/blob/main/workflows/isomiR/motif-consensus.fa)
