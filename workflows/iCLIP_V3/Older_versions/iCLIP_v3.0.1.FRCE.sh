#!/bin/bash

#SBATCH --job-name iCLIPv3.0.1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8 
#SBATCH --partition=norm
#SBATCH --time=14:00:00
#SBATCH --mem=350g

module load nextflow
module load singularity
#module load R
#module load bedtools
#module load STAR
#module load samtools
#module load umitools
#module load manorm
#module load ultraplex
#module load multiqc
#module load fastqc
#module load fastq_screen
#module load bowtie 

mkdir -p 00_QC/01_Barcodes/
mkdir -p 00_QC/02_SamStats/
mkdir -p 00_QC/03_MultiQC
mkdir -p 00_QC/04_QC_ScreenSpecies
mkdir -p 00_QC/05_QC_ScreenRRNA
mkdir -p temp
mkdir -p 01_preprocess/01_fastq/
mkdir -p 01_preprocess/02_alignment/01_unmapped/
mkdir -p 01_preprocess/07_rscripts/same
mkdir -p 01_preprocess/07_rscripts/oppo
mkdir -p 01_preprocess/07_rscripts/_STARtmp/
mkdir -p 02_bam/01_merged/
mkdir -p 02_bam/02_dedup
mkdir -p 03_peaks/01_bed/
mkdir -p 03_peaks/02_SAF/
mkdir -p 03_peaks/03_counts/
mkdir -p 04_annotation/01_project
mkdir -p 04_annotation/02_peaks/
mkdir -p 05_demethod/02_analysis
mkdir -p log/STAR

export NXF_SINGULARITY_CACHEDIR=$PWD/.singularity

timestamp=$(date +%Y%m%d_%H%M)
project="iCLIP_run_"$timestamp

Arguments=$1

nextflow run iCLIP_v3.0.1.nf -c nextflow.FRCE.config \
	--workdir $PWD \
	-params-file nextflow.parameters.FRCE.yaml \
        --outdir ${project}/ ${Arguments} \
	-with-report ${project}/Report.html \
	-with-dag ${project}/Flowchart.html \
	-with-timeline ${project}/Timeline.html


