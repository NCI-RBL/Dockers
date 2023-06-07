#!/bin/bash

#SBATCH --job-name wil_works
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=4 
#SBATCH --partition=norm,quick
#SBATCH --time=10:00:00
#SBATCH --gres=lscratch:24 
#SBATCH --mem=250g

module load nextflow
module load singularity
module load R
module load bedtools
module load STAR
module load samtools
module load umitools


# Create output file for process Check_ReadCounts:
if [[ -f !{params.workdir}/00_QC/02_SamStats/qc_read_count_raw_values.txt ]]; then rm !{params.workdir}/00_QC/02_SamStats/qc_read_count_raw_values.txt ; fi 
touch !{params.workdir}/00_QC/02_SamStats/qc_read_count_raw_values.txt

nextflow run iCLIP_light.nf -c nextflow.config --workdir $PWD --reference mm10 -params-file config/nextflow.parameters.yaml -resume
#nextflow run iCLIP_light.nf -c nextflow.config --workdir $PWD --reference mm10 -params-file config/nextflow.parameters.yaml 

