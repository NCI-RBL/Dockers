# eCLIP Workflow
The eCLIP Workflow implemented here was designed to run on a High-performance Cluster such as ***Biowulf*** or ***FRCE***. The heart of the workflow analysis uses the ***eCLIP*** workflow from Yeo Lab's.

Detailed information on required software can be found using the following links:

1. [Slurm workload Manager](https://slurm.schedmd.com)
2. [Environmental Modules](https://modules.readthedocs.io/en/latest/)
3. [SingularityCE](https://sylabs.io/singularity/)
4. [git](https://git-scm.com)
5. [eCLIP](https://github.com/YeoLab/eCLIP) from Yeo Lab.
6. [MultiQC](https://seqera.io/multiqc/) from Seqera.

The first four items are typically provided by a High-performance Cluster such as ***Biowulf*** or ***FRCE***.

## Setting Up eCLIP Workflow

### Create an Environment to Run eCLIP Workflow on FRCE

```BASH
module load mamba

mamba create -n cwl_env conda-forge::tabulate conda-forge::cwltool matplotlib

mamba activate cwl_env
pip install multiqc

multiqc --version
```

### Download Yeo Lab's GitHub Repository

```BASH
SCRIPT_DIR=/home/$USER/eCLIP_WF
mkdir $SCRIPT_DIR

cd $SCRIPT_DIR
git clone git@github.com:YeoLab/eCLIP.git

ECLIP_DIR=$SCRIPT_DIR/eCLIP
ls -l $ECLIP_DIR
```

### Setting Up Required Files and Directories

Directory locations should be based on the preferences of the user. As an example, we provide the following set up below as set up on FRCE. It is recommend that you separate the results or runs directory from the data directory. A full description of required files is located with links below:

- [Prerequisite files](https://github.com/YeoLab/eCLIP?tab=readme-ov-file#prerequisite-files)
- [Description of the manifest File](https://github.com/YeoLab/eCLIP?tab=readme-ov-file#description-of-the-manifest)

```BASH
DATA_DIR=/home/$USER/Data
RUN_DIR=/scratch/cluster_scratch/$USER/eCLIP_run

mkdir $DATA_DIR
mkdir $RUN_DIR
```

Create a manifest directory and manifest file as required by eCLIP.

```BASH
MANIFEST_DIR=$DATA_DIR/manifests
mkdir $MANIFEST_DIR

# For single-end reads copy, modify and rename file
cp $ECLIP_DIR/example/single_end_clip.yaml $MANIFEST_DIR/
```

### Modify the Run Script

Modify the run.sh script file.

## Running Workflow

The workflow has two major parts:

1. Run YeoLab's [eCLIP](https://github.com/YeoLab/eCLIP).
2. Followed by YeoLab's [merge_peaks](https://github.com/YeoLab/merge_peaks)

```BASH
sbatch run.sh
```
## Results
A full description of results files are located with [eCLIP Outputs](https://github.com/YeoLab/eCLIP?tab=readme-ov-file#outputs) and [merge_peaks Outputs](https://github.com/YeoLab/merge_peaks?tab=readme-ov-file#outputs).

We provide two summary html files are generated:

1. ***eCLIP_ReadMeSummary.html*** summary results file as a starting point for exploration of your results of ***eCLIP***, and
2. ***eCLIP_MergePeaks_ReadMeSummary.html*** for ***merge_peaks***.


