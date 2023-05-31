#!/bin/bash

#SBATCH --job-name wil_works
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4 
#SBATCH --partition=norm,quick
#SBATCH --time=12:00:00
#SBATCH --gres=lscratch:24 
#SBATCH --mem=250g

module load nextflow
module load singularity
module load R
module load bedtools

nextflow run iCLIP_light.nf -c nextflow.config --workdir $PWD --reference mm10 -params-file config/nextflow.parameters.yaml -resume

