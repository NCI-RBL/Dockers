#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --open-mode=append
#SBATCH --time=48:00:00
#SBATCH --mem=300g
#SBATCH --job-name=eCLIP
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER

module load mamba
mamba activate cwl_env

module load python
module load singularity

## Set Directory Locations
SCRIPT_DIR=/home/$USER/eCLIP_WF
ECLIP_DIR=$SCRIPT_DIR/eCLIP
RUN_DIR=/scratch/cluster_scratch/$USER/eCLIP_run/
MANIFEST_DIR=$DATA_DIR/manifests

# Switch to the Run Directory
cd $RUN_DIR

## Run pipeline
cwltool --parallel --singularity \
--tmpdir-prefix $RUN_DIR/ \
--cachedir      $RUN_DIR/ \
--leave-tmpdir \
$ECLIP_DIR/cwl/wf_get_peaks_scatter_se.cwl \
$MANIFEST_DIR/ZAP.yaml

## Run MultiQC
multiqc --title 'eCLIP' --outdir $RUN_DIR $RUN_DIR

## Run Post-processing Summary Script
clear;python $SCRIPT_DIR/src/post_process_eCLIP.py \
--yaml_file /mnt/gridftp/guibletwm/CCBRRBL16/manifests/ZAP.yaml \
--directory $RUN_DIR

## Run Merge Peaks
cwltool --parallel --singularity \
--outdir /mnt/gridftp/guibletwm/CCBRRBL16/ZAP/ \
--tmpdir-prefix /scratch/cluster_scratch/guibletwm/zap/ \
--cachedir  /scratch/cluster_scratch/guibletwm/zap/ \
--leave-tmpdir /mnt/gridftp/guibletwm/CCBRRBL16/merge_peaks/cwl/wf_get_reproducible_eclip_peaks.cwl \
/mnt/gridftp/guibletwm/CCBRRBL16/manifests/ZAP.merge.yaml

## Merge Summary
python $SCRIPT_DIR/src/summarize_merge_peaks_wf.py \
--directory $RUN_DIR \
--yaml_file /mnt/gridftp/guibletwm/CCBRRBL16/manifests/ZAP.merge.yaml


# <!DOCTYPE html>
# <html>
# <head>
#   <meta name="viewport" content="width=device-width, initial-scale=1">
#   <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/css/bootstrap.min.css">
#   <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.7.1/jquery.min.js"></script>
#   <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/js/bootstrap.min.js"></script>
# </head>
# <body>
#
# <div class="container">
#   <h2>Simple Collapsible</h2>
#   <a href="#demo" class="btn btn-info" data-toggle="collapse">Simple collapsible</a>
#   <div id="demo" class="collapse">
# 	Lorem ipsum dolor sit amet, consectetur adipisicing elit,
# 	sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam,
# 	quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.
#   </div>
# </div>
#
# </body>
# </html>
