#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=append
#SBATCH --time=2-00:00:00
#SBATCH --mem=200g
#SBATCH --job-name=DTEG
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guibletwm
  
module load nextflow
module load singularity

nextflow run -c nextflow.config test.nf --SampleInfo sample_info.txt --Process AlignRNA
nextflow run -c nextflow.config test.nf --SampleInfo sample_info.txt --Process RunRiboSeq
nextflow run -c nextflow.config test.nf --SampleInfo sample_info.txt --Process RunHTseq
nextflow run -c nextflow.config test.nf --SampleInfo sample_info.txt --Process MergeCounts
nextflow run -c nextflow.config test.nf --SampleInfo sample_info.txt --Process DTEG


