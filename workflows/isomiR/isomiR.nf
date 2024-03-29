
// Author : Wilfried Guiblet


nextflow.enable.dsl=2


// Define input arguments
params.MinLen = 18
params.fastqdir = "fastq_files/*fastq"
params.trimmer = "Qiagen"
params.outdir = "./"
params.index = "hg38"
params.consensus = ""

params.sRNAprofiling = "Yes"
params.QuagmiR = "Yes"

// Create a channel from the input path
fastq_files = Channel.fromPath(params.fastqdir)


process Trimm {
  container 'wilfriedguiblet/isomir:v0.3'

  input:
    val fastq_file

  output:
    val "${fastq_file.baseName}"

  script:
  if (params.trimmer == "Qiagen")
    """
    /opt2/adaptor_remove_Qiagen_v3 ${params.MinLen} < /data2/fastq_files/${fastq_file.baseName}.fastq > /data2/${params.outdir}/ready_files/${fastq_file.baseName}.ready.fastq 2> /data2/${params.outdir}/logs/${fastq_file.baseName}.txt
    """

  else if (params.trimmer == "Illumina")
    """
    /opt2/adaptor_remove_illumina_v5.exe ${params.MinLen} < /data2/fastq_files/${fastq_file.baseName}.fastq > /data2/${params.outdir}/ready_files/${fastq_file.baseName}.ready.fastq 2> /data2/${params.outdir}/logs/${fastq_file.baseName}.txt
    """

  else if (params.trimmer == "NEB")
    """
    /opt2/adaptor_remove_NEB_v2.exe ${params.MinLen} < /data2/fastq_files/${fastq_file.baseName}.fastq > /data2/${params.outdir}/ready_files/${fastq_file.baseName}.ready.fastq 2> /data2/${params.outdir}/logs/${fastq_file.baseName}.txt
    """

  else if (params.trimmer == "GuLab")
    """
    /opt2/adaptor_remove_lab_new8nt.exe ${params.MinLen} < /data2/fastq_files/${fastq_file.baseName}.fastq > /data2/${params.outdir}/ready_files/${fastq_file.baseName}.ready.fastq 2> /data2/${params.outdir}/logs/${fastq_file.baseName}.txt
    """

  else if (params.trimmer == "None")
    """
    cp /data2/fastq_files/${fastq_file.baseName}.fastq  /data2/${params.outdir}/ready_files/${fastq_file.baseName}.ready.fastq
    """

}



process sRNAprofiling {

  container 'wilfriedguiblet/isomir:v0.3'

  input:
    val fastq_ready

  output:
    val fastq_ready

  script:

  if (params.sRNAprofiling == "Yes")
  """
  echo 'processing ${fastq_ready}..'
  echo 'mapping to rRNA..mode -n 0 -l 20...'
  bowtie -n 0 -l 20 --best --norc /data2/index/${params.index}/rRNA --un /data2/${params.outdir}/analysis_files/'${fastq_ready}'_rm_rRNA.fastq --al /data2/${params.outdir}/analysis_files/'${fastq_ready}'_rRNA.fastq /data2/${params.outdir}/ready_files/'${fastq_ready}'.ready.fastq /data2/${params.outdir}/analysis_files/temp.txt

  echo 'mapping to tRNA..mode -n 0 -l 20...'
  bowtie -n 0 -l 20 --best --norc /data2/index/${params.index}/tRNA --un /data2/${params.outdir}/analysis_files/'${fastq_ready}'_rm_rRNA_tRNA.fastq --al /data2/${params.outdir}/analysis_files/'${fastq_ready}'_tRNA.fastq /data2/${params.outdir}/analysis_files/'${fastq_ready}'_rm_rRNA.fastq /data2/${params.outdir}/analysis_files/temp.txt

  echo 'mapping to snoRNA..mode -n 0 -l 20...'
  bowtie -n 0 -l 20 --best --norc /data2/index/${params.index}/snoRNA --un /data2/${params.outdir}/analysis_files/'${fastq_ready}'_rm_rRNA_tRNA_snoRNA.fastq --al /data2/${params.outdir}/analysis_files/'${fastq_ready}'_snoRNA.fastq /data2/${params.outdir}/analysis_files/'${fastq_ready}'_rm_rRNA_tRNA.fastq /data2/${params.outdir}/analysis_files/temp.txt

  echo 'mapping to miRNA_precusor.. n mode mismatch 0 seed length 20...'
  bowtie -n 0 -l 20 --best --norc --sam /data2/index/${params.index}/miRNA_hairpin_all --un /data2/${params.outdir}/analysis_files/'${fastq_ready}'_rm_rRNA_tRNA_snoRNA_miRNA.fastq --al /data2/${params.outdir}/analysis_files/'${fastq_ready}'_miRNA.fastq /data2/${params.outdir}/analysis_files/'${fastq_ready}'_rm_rRNA_tRNA_snoRNA.fastq /data2/${params.outdir}/analysis_files/'${fastq_ready}'_miRNA.sam
  bowtie -n 0 -l 20 --best --norc --sam /data2/index/${params.index}/miRNA_hairpin_hc --al /data2/${params.outdir}/analysis_files/'${fastq_ready}'_hc_miRNA.fastq /data2/${params.outdir}/analysis_files/'${fastq_ready}'_miRNA.fastq /data2/${params.outdir}/analysis_files/'${fastq_ready}'_hc_miRNA.sam

  echo 'mapping to mRNA.. mode -n 0 -l 20...'
  bowtie -n 0 -l 20 --best --norc /data2/index/${params.index}/refMrna --un /data2/${params.outdir}/analysis_files/'${fastq_ready}'_rm_rRNA_tRNA_snoRNA_miRNA_mRNA.fastq --al /data2/${params.outdir}/analysis_files/'${fastq_ready}'_mRNA.fastq /data2/${params.outdir}/analysis_files/'${fastq_ready}'_rm_rRNA_tRNA_snoRNA_miRNA.fastq /data2/${params.outdir}/analysis_files/temp.txt

  echo 'mapping to refSeq.. mode -n 0 -l 20...'
  bowtie -n 0 -l 20 --best --norc /data2/index/${params.index}/ref_transcripts --un /data2/${params.outdir}/analysis_files/'${fastq_ready}'_others.fastq --al /data2/${params.outdir}/analysis_files/'${fastq_ready}'_others_ref.fastq /data2/${params.outdir}/analysis_files/'${fastq_ready}'_rm_rRNA_tRNA_snoRNA_miRNA_mRNA.fastq /data2/${params.outdir}/analysis_files/temp.txt

  echo 'processing ${fastq_ready}...'
  echo 'mapping to mycoplasma H..mode v1...'
  bowtie -v 1 --sam /data2/index/mycoplasma/Mycoplasma_hyo --un /data2/${params.outdir}/analysis_files/'${fastq_ready}'_unmappable.fastq  --al /data2/${params.outdir}/analysis_files/'${fastq_ready}'_mycoplasmaH.fastq /data2/${params.outdir}/analysis_files/'${fastq_ready}'_others.fastq  /data2/${params.outdir}/analysis_files/'${fastq_ready}'_mycoplasma.sam

  total=`cat /data2/${params.outdir}/ready_files/${fastq_ready}.ready.fastq | wc -l`
  rRNA=`cat /data2/${params.outdir}/analysis_files/${fastq_ready}_rRNA.fastq | wc -l`
  tRNA=`cat /data2/${params.outdir}/analysis_files/${fastq_ready}_tRNA.fastq | wc -l`
  snoRNA=`cat /data2/${params.outdir}/analysis_files/${fastq_ready}_snoRNA.fastq | wc -l`
  miRNA=`cat /data2/${params.outdir}/analysis_files/${fastq_ready}_miRNA.fastq | wc -l`
  mRNA=`cat /data2/${params.outdir}/analysis_files/${fastq_ready}_mRNA.fastq | wc -l`
  others_ref=`cat /data2/${params.outdir}/analysis_files/${fastq_ready}_others_ref.fastq | wc -l`
  mycoplasma=`cat /data2/${params.outdir}/analysis_files/${fastq_ready}_mycoplasmaH.fastq | wc -l`
  unmappable=`cat /data2/${params.outdir}/analysis_files/${fastq_ready}_unmappable.fastq | wc -l`
  hc_miRNA=`cat /data2/${params.outdir}/analysis_files/${fastq_ready}_hc_miRNA.fastq | wc -l`

  let total=\${total}/4
  let rRNA=\${rRNA}/4
  let tRNA=\${tRNA}/4
  let snoRNA=\${snoRNA}/4
  let miRNA=\${miRNA}/4
  let mRNA=\${mRNA}/4
  let others_ref=\${others_ref}/4
  let mycoplasma=\${mycoplasma}/4
  let unmappable=\${unmappable}/4
  let hc_miRNA=\${hc_miRNA}/4

  printf "%s\t%d\t%d\t%d\t%d\t%d\t%d     \t%d\t%d\t%d\t%d\n" ${fastq_ready} \$total \$rRNA \$tRNA \$snoRNA \$miRNA \$mRNA \$others_ref \$mycoplasma \$unmappable \$hc_miRNA >> /data2/${params.outdir}/analysis_results/small_RNA_profile.txt
  """

}


process RunQuagmiR {

  container 'wilfriedguiblet/isomir:v0.3'

  input:
    val fastq_ready

  output:
    val fastq_ready

  script:



  if(params.consensus == '')
    """
    cp -r /opt2/QuagmiR/ /data2/${params.outdir}/
    rm /data2/${params.outdir}/QuagmiR/data/sample.fastq
    ln -s /data2/${params.outdir}/ready_files/*.ready.fastq /data2/${params.outdir}/QuagmiR/data/
    cd /data2/${params.outdir}/QuagmiR/
    snakemake -j
    """

  else if(params.consensus != '')
    """
    cp -r /opt2/QuagmiR/ /data2/${params.outdir}/
    rm /data2/${params.outdir}/QuagmiR/data/sample.fastq
    cp /data2/${params.consensus} /data2/${params.outdir}/QuagmiR/motif-consensus.fa
    ln -s /data2/${params.outdir}/ready_files/*.ready.fastq /data2/${params.outdir}/QuagmiR/data/
    cd /data2/${params.outdir}/QuagmiR/
    snakemake -j
    """


}

process RPM_summary {

  container 'rocker/tidyverse'

  input:
    val fastq_ready

  script:
    """
    Rscript /data2/${params.outdir}/QuagmiR/RPM_summary.R /data2/${params.outdir}/QuagmiR/group_results/cohort1.isomir.tsv /data2/${params.outdir}/RPM_summary.csv
    """
}


workflow {

  if (params.sRNAprofiling == "Yes" && params.QuagmiR == "Yes") {Trimm(fastq_files) | sRNAprofiling | collect | RunQuagmiR | RPM_summary}
  if (params.sRNAprofiling == "No" && params.QuagmiR == "Yes") {Trimm(fastq_files)  | collect | RunQuagmiR | RPM_summary}
  if (params.sRNAprofiling == "Yes" && params.QuagmiR == "No") {Trimm(fastq_files) | sRNAprofiling}
  if (params.sRNAprofiling == "No" && params.QuagmiR == "No") {Trimm(fastq_files)}

}














