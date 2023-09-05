#!/bin/bash
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --open-mode=append
#SBATCH --time=1:00:00
#SBATCH --mem=150g
#SBATCH --job-name=RiboPrint

# Author : Wilfried Guiblet

module load singularity
module load nextflow
module load ucsc
module load bedtools
module load R

export NXF_SINGULARITY_CACHEDIR=$PWD/.singularity
export SINGULARITY_CACHEDIR=$PWD/.singularity

mkdir -p Results

ResumeArg=$1

nextflow run RiboFootPrint.nf --workdir $PWD  -c nextflow.config -params-file RiboFootPrint.parameters.yaml  ${ResumeArg} 
#nextflow run RiboFootPrint.nf -c nextflow.config --outdir Results/ --workdir $PWD -resume 

#gzip Results/*bedgraph


