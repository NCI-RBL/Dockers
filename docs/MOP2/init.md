### Prepare directories

Create a working directory and subdirectories:

```bash
mkdir MOP2
cd MOP2
mkdir MOP2_work
mkdir MOP2_output
```

### Scripts and dependencies

Copy the following files and folders in your working directory:

```bash
cp /scratch/cluster_scratch/guibletwm/MOP2_repo/ ./
```


Download the following files in the working directory:

```bash
- [PolyATail.sh](https://github.com/RBL-NCI/Dockers/blob/main/workflows/MOP2/PolyATail.sh)
```


### Paths

Modify the following paths in PolyATail.sh:

- basedir
- nextflow
- procdir
