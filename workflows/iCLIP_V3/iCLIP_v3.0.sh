#!/bin/bash

#SBATCH --job-name wil_works
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2 
#SBATCH --partition=norm,quick
#SBATCH --time=8:00:00
#SBATCH --gres=lscratch:24 
#SBATCH --mem=250g

module load nextflow
module load singularity
module load R
module load bedtools
module load STAR
module load samtools
module load umitools
module load manorm
module load ultraplex
module load multiqc
module load fastqc
module load fastq_screen
module load bowtie 

nextflow run iCLIP_v3.0.nf -c nextflow.config \
	--workdir $PWD \
	--reference mm10 \
	-params-file config/nextflow.parameters.yaml \
        --outdir ${project}/ ${Arguments} \
	-with-report ${project}/Report.html \
	-with-dag ${project}/Flowchart.html \
	-with-timeline ${project}/Timeline.html \
	-resume
