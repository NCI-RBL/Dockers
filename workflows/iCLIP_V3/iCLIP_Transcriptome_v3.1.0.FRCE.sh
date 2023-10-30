#!/bin/bash

#SBATCH --job-name iCLIPv3.1.0
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8 
#SBATCH --partition=ccr,norm
#SBATCH --time=18:00:00
#SBATCH --mem=350g

module load nextflow
module load singularity


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
export SINGULARITY_CACHEDIR=$PWD/.singularity

timestamp=$(date +%Y%m%d_%H%M)
project="iCLIP_run_"$timestamp

Arguments=$1


nextflow run iCLIP_Transcriptome_v3.1.0.nf -c nextflow.FRCE.config \
        --workdir $PWD \
        -params-file nextflow.parameters.FRCE.yaml \
        -with-report ${project}/Report.html \
        -with-dag ${project}/Flowchart.html \
        -with-timeline ${project}/Timeline.html \
        ${Arguments}



