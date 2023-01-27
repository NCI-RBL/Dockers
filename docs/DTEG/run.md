Simply submit the pipeline as a slurm job:

sbatch DTEG_RBL.sh


If you do not have slurm and/or are running the pipeline on a local computer (not recommended), you may run the different parts of the workflow as follow:

nextflow run -c nextflow.config test.nf --SampleInfo sample_info.txt --Process AlignRNA
nextflow run -c nextflow.config test.nf --SampleInfo sample_info.txt --Process RunRiboSeq
nextflow run -c nextflow.config test.nf --SampleInfo sample_info.txt --Process RunHTseq
nextflow run -c nextflow.config test.nf --SampleInfo sample_info.txt --Process MergeCounts
nextflow run -c nextflow.config test.nf --SampleInfo sample_info.txt --Process DTEG


Running this pipeline will download the Docker container in a folder named "work". To remove the container, delete the folder.

rm -r work
