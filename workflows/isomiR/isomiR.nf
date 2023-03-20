
// Author : Wilfried Guiblet


nextflow.enable.dsl=2


// Define input arguments
params.MinLen = 18
params.fastqdir = "fastq_files/*fastq"


// Create a channel from the input path
fastq_files = Channel.fromPath(params.fastqdir)


process Trimm {
  container 'wilfriedguiblet/isomir:v0.2'

  input:
    val fastq_file

  output:
    val "${fastq_file.baseName}"

  script:
  """
  /opt2/adaptor_remove_Qiagen_v3 ${params.MinLen} < /data2/fastq_files/${fastq_file.baseName}.fastq > /data2/ready_files/${fastq_file.baseName}.ready.fastq
  """
}



process sRNAprofiling {

  container 'wilfriedguiblet/isomir:v0.2'

  input:
    val fastq_ready

  output:
    val fastq_ready

  script:
  """
  echo 'processing ${fastq_ready}..'
  echo 'mapping to rRNA..mode -n 0 -l 20...'
  bowtie -n 0 -l 20 --best --norc /data2/NGS/index/human/rRNA --un /data2/analysis_files/'${fastq_ready}'_rm_rRNA.fastq --al /data2/analysis_files/'${fastq_ready}'_rRNA.fastq /data2/ready_files/'${fastq_ready}'.ready.fastq /data2/analysis_files/temp.txt

  echo 'mapping to tRNA..mode -n 0 -l 20...'
  bowtie -n 0 -l 20 --best --norc /data2/NGS/index/human/tRNA --un /data2/analysis_files/'${fastq_ready}'_rm_rRNA_tRNA.fastq --al /data2/analysis_files/'${fastq_ready}'_tRNA.fastq /data2/analysis_files/'${fastq_ready}'_rm_rRNA.fastq /data2/analysis_files/temp.txt

  echo 'mapping to snoRNA..mode -n 0 -l 20...'
  bowtie -n 0 -l 20 --best --norc /data2/NGS/index/human/snoRNA --un /data2/analysis_files/'${fastq_ready}'_rm_rRNA_tRNA_snoRNA.fastq --al /data2/analysis_files/'${fastq_ready}'_snoRNA.fastq /data2/analysis_files/'${fastq_ready}'_rm_rRNA_tRNA.fastq /data2/analysis_files/temp.txt

  echo 'mapping to miRNA_precusor.. n mode mismatch 0 seed length 20...'
  bowtie -n 0 -l 20 --best --norc --sam /data2/NGS/index/miRNA_Aug_2018/hsa_hairpin_all_2018 --un /data2/analysis_files/'${fastq_ready}'_rm_rRNA_tRNA_snoRNA_miRNA.fastq --al /data2/analysis_files/'${fastq_ready}'_miRNA.fastq /data2/analysis_files/'${fastq_ready}'_rm_rRNA_tRNA_snoRNA.fastq /data2/analysis_files/'${fastq_ready}'_miRNA.sam
  bowtie -n 0 -l 20 --best --norc --sam /data2/NGS/index/miRNA_Aug_2018/hsa_hairpin_hc_2018 --al /data2/analysis_files/'${fastq_ready}'_hc_miRNA.fastq /data2/analysis_files/'${fastq_ready}'_miRNA.fastq /data2/analysis_files/'${fastq_ready}'_hc_miRNA.sam

  echo 'mapping to mRNA.. mode -n 0 -l 20...'
  bowtie -n 0 -l 20 --best --norc /data2/NGS/index/human/refMrna --un /data2/analysis_files/'${fastq_ready}'_rm_rRNA_tRNA_snoRNA_miRNA_mRNA.fastq --al /data2/analysis_files/'${fastq_ready}'_mRNA.fastq /data2/analysis_files/'${fastq_ready}'_rm_rRNA_tRNA_snoRNA_miRNA.fastq /data2/analysis_files/temp.txt

  echo 'mapping to refSeq.. mode -n 0 -l 20...'
  bowtie -n 0 -l 20 --best --norc /data2/NGS/index/human/ref_transcripts --un /data2/analysis_files/'${fastq_ready}'_others.fastq --al /data2/analysis_files/'${fastq_ready}'_others_ref.fastq /data2/analysis_files/'${fastq_ready}'_rm_rRNA_tRNA_snoRNA_miRNA_mRNA.fastq /data2/analysis_files/temp.txt

  echo 'processing ${fastq_ready}...'
  echo 'mapping to mycoplasma H..mode v1...'
  bowtie -v 1 --sam /data2/NGS/index/Mycoplasma_hyo --un /data2/analysis_files/'${fastq_ready}'_unmappable.fastq  --al /data2/analysis_files/'${fastq_ready}'_mycoplasmaH.fastq /data2/analysis_files/'${fastq_ready}'_others.fastq  /data2/analysis_files/'${fastq_ready}'_mycoplasma.sam

  total=`cat /data2/ready_files/${fastq_ready}.ready.fastq | wc -l`
  rRNA=`cat /data2/analysis_files/${fastq_ready}_rRNA.fastq | wc -l`
  tRNA=`cat /data2/analysis_files/${fastq_ready}_tRNA.fastq | wc -l`
  snoRNA=`cat /data2/analysis_files/${fastq_ready}_snoRNA.fastq | wc -l`
  miRNA=`cat /data2/analysis_files/${fastq_ready}_miRNA.fastq | wc -l`
  mRNA=`cat /data2/analysis_files/${fastq_ready}_mRNA.fastq | wc -l`
  others_ref=`cat /data2/analysis_files/${fastq_ready}_others_ref.fastq | wc -l`
  mycoplasma=`cat /data2/analysis_files/${fastq_ready}_mycoplasmaH.fastq | wc -l`
  unmappable=`cat /data2/analysis_files/${fastq_ready}_unmappable.fastq | wc -l`
  hc_miRNA=`cat /data2/analysis_files/${fastq_ready}_hc_miRNA.fastq | wc -l`

  let total=\${total}/4
  let rRNA=\${rRNA}/4
  let tRNA=\${tRNA}/4
  let snoRNA=\${snoRNA}/4
  let miRNA=\${miRNA}/4
  let mRNA=\${mRNA}/4
  let others_ref=\${others_ref}/4
  let mycoplasma=\${mycoplasma}/4
  let unmappable=\${unmappable}/4

  printf "%s\t%d\t%d\t%d\t%d\t%d\t%d     \t%d\t%d\t%d\t%d\n" ${fastq_ready} \$total \$rRNA \$tRNA \$snoRNA \$miRNA \$mRNA \$others_ref \$mycoplasma \$unmappable \$hc_miRNA >> /data2/analysis_results/small_RNA_profile.txt
  """

}


process RunQuagmiR {

  container 'wilfriedguiblet/isomir:v0.2'

  input:
    val fastq_ready

  script:
  """
  cp -r /opt2/QuagmiR/ /data2/
  ln -s /data2/ready_files/${fastq_ready}.ready.fastq /data2/QuagmiR/data/
  cd /data2/QuagmiR/
  snakemake -j    
  """
}




workflow {
  Trimm(fastq_files) | sRNAprofiling | RunQuagmiR


}













