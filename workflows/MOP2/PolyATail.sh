#!/bin/bash
 
# Location of software installation
 
basedir="/mnt/rnabl-work/Guiblet/CCBRRBL7/MOP2_repo"
nextflow="/mnt/rnabl-work/Guiblet/CCBRRBL7/nextflow"
 
# Location of reference genome and annotations
 
ref="/mnt/RBL-MRRD-CryoEM/static/references_genomes/gencode.v38.transcripts.fa"
annot="/mnt/RBL-MRRD-CryoEM/static/references_genomes/gencode.v38.annotation.gtf"

# location of processing directory
 
procdir="/mnt/rnabl-work/Guiblet/CCBRRBL7/MOP2_work/"

# Define location of Nanopore sequencing data in Fast5 format
 
echo -e "Using data found in $1"
 
fast5="$1"
    
# examining passed/failed reads in the specified location
 
fast5_total_num=`find $fast5 -type f -name '*.fast5' | sort -u | awk -F "/" '{ print $0}' | wc -l`
fast5_pass_num=`find $fast5 -type f -name '*.fast5' | sort -u | awk -F "/" '{ print $0}' | awk '/(pass)/' | wc -l`
fast5_fail_num=`find $fast5 -type f -name '*.fast5' | sort -u | awk -F "/" '{ print $0}' | awk '/(fail)/' | wc -l`
 
echo -e "Found $fast5_total_num reads with $fast5_pass_num passed and $fast5_fail_num failed"

echo -ne "Provide a short project name: "
read project_name

# create processing directories
 
timestamp=$(date +%Y%m%d_%H%M)
 
project=""$procdir""$timestamp"_master_of_pores_V2_"$project_name""
 
mkdir $project
mkdir $project/data

# create a folder with soft links to the reads
 
mkdir $project/data/fast5_dir
 
for file in `find $fast5 -type f -name '*.fast5' | sort -u | awk -F "/" '{ print $0}' | awk '!/(skip|fail)/'`; do
                read_name=`echo $file | sed 's!.*/!!'`;
                ln -s $file $project/data/fast5_dir/$read_name;
done
 
  
# location of output
 
outdir="/mnt/rnabl-work/Guiblet/CCBRRBL7/MOP2_output/"$timestamp"_master_of_pores_V2_"$project_name""
mkdir ${outdir}
mkdir ${outdir}/mop_preprocess
mkdir ${outdir}/mop_tail
  
# making further processing directories and hard links
 
ln $basedir/nextflow.global.config $project/nextflow.global.config
mkdir $project/singularity
ln $basedir/singularity/* $project/singularity/
mkdir $project/mop_preprocess/
ln $basedir/mop_preprocess/mop_preprocess.nf $project/mop_preprocess/mop_preprocess.nf
ln $basedir/mop_preprocess/config.yaml $project/mop_preprocess/config.yaml
ln $basedir/mop_preprocess/nextflow.config $project/mop_preprocess/nextflow.config
mkdir $project/mop_preprocess/deeplexicon
ln $basedir/mop_preprocess/deeplexicon/* $project/mop_preprocess/deeplexicon/
mkdir $project/mop_preprocess/bin
cp -r $basedir/mop_preprocess/bin/ont-guppy/ $project/mop_preprocess/bin/
ln -s $project/mop_preprocess/bin/ont-guppy/bin/*guppy* $project/mop_preprocess/bin/
ln $basedir/mop_preprocess/bin/*.py $project/mop_preprocess/bin/
mkdir $project/conf
ln $basedir/conf/* $project/conf/
ln $basedir/outdirs.nf $project/
ln $basedir/local_modules.nf $project/
ln $basedir/mop_preprocess/drna_tool_splice_opt.tsv $project/mop_preprocess/drna_tool_splice_opt.tsv
ln $basedir/mop_preprocess/final_summary_01.txt $project/mop_preprocess/ # DOES THIS FILE NEED UPDATE WITH EACH PROJECT????
mkdir $project/docs/
ln $basedir/docs/* $project/docs/
mkdir $project/img
ln $basedir/img/* $project/img/
cp -r $basedir/BioNextflow/ $project/

echo -e "Directories created and files copied"
 
### Pre-processing ###
 
echo "\
params {

    kit                 = \"SQK-RNA002\"
    flowcell            = \"FLO-MIN106\"

    conffile            = \"final_summary_01.txt\"
    fast5               = \"${project}/data/fast5_dir/*.fast5\"
    fastq               = \"\"

    reference           = \"$ref\"
    annotation          = \"$annot\"
    ref_type            = \"transcriptome\"

    pars_tools          = \"drna_tool_splice_opt.tsv\"
    seq_type            = \"RNA\"
    output              = \"${project}/mop_preprocess/output/\"
    qualityqc           = 5
    granularity         = 1

    basecalling          = \"guppy\"
    basecalling_opt      = \"\"
    GPU                 = \"ON\"
    demultiplexing      = \"NO\"
    demultiplexing_opt  = \"NO\"
    demulti_fast5       = \"NO\"

    filtering              = \"nanoq\"
    filtering_opt          = \"\"

    mapping              = \"graphmap\"
    mapping_opt          = \"minimap2\"
    map_type            = \"spliced\"

    counting             = \"nanocount\"
    counting_opt         = \"\"

    variant_calling      = \"NO\"
    variant_opt         = \"\"    

    discovery           = \"NO\"

    cram_conv           = \"YES\"
    subsampling_cram    = 50

    saveSpace           = \"NO\"

 
    downsampling        = \"\"
 
    email               = \"\"
}" > ${project}/mop_preprocess/params.config
 
echo "\
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=32
#SBATCH --open-mode=append
#SBATCH --time=2-00:00:00
#SBATCH --mem=20g
#SBATCH --job-name=MOP2_PREPROCESS
#SBATCH --output=${project}/mop_preprocess/mop_preprocess.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guibletwm
  
hostname
nvidia-smi
module purge
module load singularity/3.7.2
module load nextflow 

cd $project/mop_preprocess
nextflow run mop_preprocess.nf -with-singularity
ln -s "$outdir"/mop_preprocess RNA
cp params.config "$outdir"/mop_preprocess/
cp mop_preprocess.log "$outdir"/mop_preprocess/
cp masterofpores_preprocess.sh "$outdir"/mop_preprocess/
" > $project/mop_preprocess/masterofpores_preprocess.sh
 
sbatch ${project}/mop_preprocess/masterofpores_preprocess.sh
 
### PolyA Tail analysis ####

# making further processing directories and hard links
 
mkdir $project/mop_tail
 
ln $basedir/mop_tail/mop_tail.nf $project/mop_tail/mop_tail.nf
ln $basedir/mop_tail/tools_opt.tsv $project/mop_tail/tools_opt.tsv
ln $basedir/mop_tail/nextflow.config $project/mop_tail/nextflow.config
mkdir $project/mop_tail/bin
ln $basedir/mop_tail/bin/* $project/mop_tail/bin/

echo -e "Directories created and files copied"
 
### Computing PolyATails ###

echo "\
params {
    input_path      = \"../mop_preprocess/output/\"
    reference          = \"$ref\"
    pars_tools         = \"tools_opt.tsv\"

    tailfindr          = \"YES\"
    nanopolish         = \"YES\"

    nanopolish_opt     = \"\"
    tailfindr_opt      = \"basecall_group = 'Basecall_1D_000'\"
 
    output             = \"$project/mop_tail/outputPoly/\"
 
    email              = \"\"
}" > $project/mop_tail/params.config
 
 
echo "\
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=append
#SBATCH --time=2-00:00:00
#SBATCH --mem=20g
#SBATCH --job-name=MOP2_TAIL
#SBATCH --output=$project/mop_tail/mop_tail.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guibletwm
 
hostname
module purge
module load R
module load singularity/3.7.2
module load nextflow
  
cd $project/mop_tail
 
while [ 1 ]; do
    [ ! -e ../mop_preprocess/output/report/multiqc_report.html ] && sleep 1 || break
done
 
nextflow run mop_tail.nf -with-singularity -profile slurm
cp params.config "$outdir"/mop_tail/
cp mop_tail.log "$outdir"/mop_tail/
cp masterofpores_tail.sh "$outdir"/mop_tail/
cp -r outputPoly "$outdir"/mop_tail/
" > $project/mop_tail/masterofpores_tail.sh
 
sbatch $project/mop_tail/masterofpores_tail.sh




