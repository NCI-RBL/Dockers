#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=append
#SBATCH --time=12:00:00
#SBATCH --mem=200g
#SBATCH --job-name=RiboPrint


# Author : Wilfried Guiblet

module load singularity
module load nextflow

nextflow run RiboFootPrint.nf -c nextflow.config --outdir Results/ -resume 
