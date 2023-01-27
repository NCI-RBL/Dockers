
// Author : Wilfried Guiblet


// Define empty lists
def Ribo = []
def RNA = []


// Parse the SampleInfo file and extract the Ribo and RNA Sample IDs
file(params.SampleInfo)
    .readLines()
    .each { def fields = it.split("\t")
           if (fields[2].contains('RIBO')) {Ribo.add(fields[0])}
           if (fields[2].contains('RNA')) {RNA.add(fields[0])}
          }


//Ribo = ['HEK_WT1','HEK_WT2','HEK_ThumD_KO1','HEK_ThumD_KO2']
//RNA = ['HEK293WT1RNA','HEK293WT2RNA','HEK293ThumpD1KO1RNA','HEK293ThumpD1KO2RNA']


nextflow.enable.dsl=2

params.Ribo_join = Ribo.join("_counts.txt ") + "_counts.txt"
params.RNA_join = RNA.join("_counts.txt ") + "_counts.txt"


process AlignRNA {
  container 'wilfriedguiblet/dteg:v0.3'

  input:
    val infile 
  output:
    stdout
  """
  StarAlign_RNAseq_paired.py /data2/FASTQ/ ${infile} /data2/index/hg19/ 1
  """
}



params.datapath = "/data2/FASTQ/"
params.fastq_filein = Ribo.join(",")
params.exp_name = Ribo.join(",")
params.ncrnaDir = "/data2/index/star_hg19_ncrna/"
params.genomeDir = "/data2/index/hg19/"
params.GTFfile = "/data2/index/hg19/"
params.twobitfile = "/data2/index/hg19/hg19.2bit"
params.UTRfilestring = "/data2/index/hg19/CanonicalTranscriptKnownGeneCoding_UTRs.csv"
params.mismatch = "0.1"
params.fiveprimetrim = "4"
params.linker = "NNNNNNCACTCGGGCACCAAGGA"
params.threshold = "0"
params.filtermodule = "0"
params.exclustionmodule = "0"
params.equalweight = "1"
params.ftsize = "'15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45'"
params.regionlength5 = "50"
params.regionlength3 = "100"
params.normalization = "1"
params.Rpath = "/data2/R/"
params.sizedistr_outpath = "/data2/size_distr"
params.outpath = "/data2/output/"


process RunRiboSeq {
  output:
    stdout
  container 'wilfriedguiblet/dteg:v0.3'
  """
  #80S_basic_workflow.py /data2/settings.txt
  80S_basic_workflow.py --datapath ${params.datapath} \
                        --fastq_filein ${params.fastq_filein} \
                        --exp_name ${params.exp_name} \
                        --ncrnaDir ${params.ncrnaDir} \
                        --genomeDir ${params.genomeDir} \
                        --GTFfile ${params.GTFfile} \
                        --twobitfile ${params.twobitfile} \
                        --UTRfilestring ${params.UTRfilestring} \
                        --mismatch ${params.mismatch} \
                        --fiveprimetrim ${params.fiveprimetrim} \
                        --linker ${params.linker} \
                        --threshold ${params.threshold} \
                        --filtermodule ${params.filtermodule} \
                        --exclustionmodule ${params.exclustionmodule} \
                        --equalweight ${params.equalweight} \
                        --ftsize ${params.ftsize} \
                        --regionlength5 ${params.regionlength5} \
                        --regionlength3 ${params.regionlength3} \
                        --normalization ${params.normalization} \
                        --Rpath ${params.Rpath} \
                        --sizedistr_outpath ${params.sizedistr_outpath} \
                        --outpath ${params.outpath}

  """
}

process RunHTseq {
  input:
    val infile
  output:
    stdout
  container 'wilfriedguiblet/dteg:v0.3'
  """
  HTSeq_RNAseq.py /data2/FASTQ/${infile}/${infile}_star*/${infile}_match.bam /data2/index/hg19/CanonicalTranscriptKnownGeneCoding.gtf /data2/HTSeq/${infile}_counts.txt
  """
}


process MergeCounts {
  output:
    stdout
  container 'wilfriedguiblet/dteg:v0.3'
  """
  Merge_HTSeq_Counts.py /data2/HTSeq/ RiboSeqCounts.txt ${params.Ribo_join}
  Merge_HTSeq_Counts.py /data2/HTSeq/ RNASeqCounts.txt ${params.RNA_join}
  """

}


process DTEG{
  output:
    stdout
  container 'wilfriedguiblet/dteg:v0.3'
  """
  DTEG.R /data2/HTSeq/RiboSeqCounts.txt /data2/HTSeq/RNASeqCounts.txt /data2/sample_info.txt 1
  """
}

workflow {
  if (params.Process.contains('AlignRNA')) {Channel.fromList(RNA) | AlignRNA | view}
  if (params.Process.contains('RunRiboSeq')) {RunRiboSeq | view}
  if (params.Process.contains('RunHTseq')) {Channel.fromList(Ribo + RNA) | RunHTseq | view}
  if (params.Process.contains('MergeCounts')) {MergeCounts | view}
  if (params.Process.contains('DTEG')) {DTEG | view}
}













