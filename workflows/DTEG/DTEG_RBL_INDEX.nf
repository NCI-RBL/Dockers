process genomeGenerate {
  container 'wilfriedguiblet/dteg:v0.3'
  '''
  STAR --runThreadN 20 --runMode genomeGenerate --genomeSAindexNbases 7 --genomeSAsparseD 1 --genomeDir /data2/index/star_hg19_ncrna --genomeFastaFiles /data2/index/star_hg19_ncrna/hg19_ncrna.fa --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon CDS
  STAR --runThreadN 20 --runMode genomeGenerate --genomeDir /data2/index/hg19/ --genomeFastaFiles /data2/index/hg19/hg19.fa --sjdbGTFfile /data2/index/hg19/KnownGene.gtf
  faToTwoBit /data2/index/hg19_PPL/hg19_PPL.fa /data2/index/hg19_PPL/hg19_PPL.2bit

  '''
}

workflow {
   genomeGenerate | view { it.trim() }
}

