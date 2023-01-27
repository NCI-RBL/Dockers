### Description

This document describes how to use the Master Of Pores V2 pipeline on the FRCE server.


### Initialize

Copy the following files and folders in your directory:

/scratch/cluster_scratch/guibletwm/MOP2_repo/
/scratch/cluster_scratch/guibletwm/nextflow
/scratch/cluster_scratch/guibletwm/PolyATail.sh

Create the following empty folders:

MOP2_work
MOP2_output


Modify the following paths in PolyATail.sh:

basedir
nextflow
procdir


### Run pipeline

Use the command:

bash PolyATail.sh /path_to_fast5_files/

The script will ask for a short project name. Intermediary files will be stored in the MOP2_work folder. PolyA-tail lengths will be copied in the MOP2_output.
