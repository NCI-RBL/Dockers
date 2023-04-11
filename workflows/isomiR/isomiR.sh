#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=append
#SBATCH --time=12:00:00
#SBATCH --mem=200g
#SBATCH --job-name=isomiR
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guibletwm


# Author : Wilfried Guiblet

set -eu

project_name=$1
trimmer=$2
MinLen=$3
Consensus=$4
procdir=./

timestamp=$(date +%Y%m%d_%H%M)
 
project=""$procdir""$timestamp"_isomiR_"$project_name""
 
mkdir $project
mkdir ${project}/ready_files
mkdir ${project}/analysis_files
mkdir ${project}/analysis_results


#path=`pwd`/${project}

#echo "\
#singularity {
#        enabled = true
#        autoMounts = true
#        runOptions = '"-B ${path}:/data2/"'
#}" > ./nextflow.config

#echo "\
#singularity {
#        enabled = true
#        autoMounts = true
#        runOptions = "-B $PWD:/data2/"
#}" > ./nextflow.config


printf "\ttotal\trRNA\ttRNA\tsnoRNA\tmiRNA\tmRNA\tothers_ref\tmycoplasma_H\tunmappable\thc_miRNA\n" > ${project}/analysis_results/small_RNA_profile.txt

module load nextflow
module load singularity

nextflow run -c nextflow.config isomiR.nf --outdir ${project}/ --trimmer ${trimmer} --MinLen ${MinLen} --consensus ${Consensus}

