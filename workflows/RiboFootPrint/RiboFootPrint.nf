// Author : Wilfried Guiblet

nextflow.enable.dsl=2


// Define input arguments
params.fastqdir = 'Project/40S/'
params.fastq_files = 'Project/40S/*fastq.gz'
params.outdir = ''
params.linker = 'NNNNNNCACTCGGGCACCAAGGA'

// Create a channel from the input path
fastq_files = Channel.fromPath(params.fastq_files)


process RunRiboSeq {

  container 'wilfriedguiblet/dteg:v0.3'

  input:
    val fastq_file

  output:
    val "${fastq_file.baseName}"


  script:
    """
    seqtk trimfq -b 4 /data2/${params.fastqdir}/${fastq_file.baseName}.gz  > /data2/${params.outdir}/${fastq_file.baseName}.trimmed
    skewer -t 26 -x ${params.linker} -l 15 /data2/${params.outdir}/${fastq_file.baseName}.trimmed

    STAR --runThreadN 26 --readFilesIn /data2/${params.outdir}/${fastq_file.baseName}-trimmed.fastq --genomeDir /data2/index/star_hg19_ncrna/ --outReadsUnmapped Fastx --outSAMtype None --outFilterMismatchNoverLmax 0.3 --outFileNamePrefix /data2/${params.outdir}/${fastq_file.baseName}.ncrnaout.

    STAR --runThreadN 26 --readFilesIn /data2/${params.outdir}/${fastq_file.baseName}.ncrnaout.Unmapped.out.mate1 --limitBAMsortRAM 600000000000 --genomeDir /data2/index/hg19/ --outFilterType BySJout --outFilterIntronMotifs RemoveNoncanonicalUnannotated --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.1 --outWigType wiggle read1_5p --outWigNorm RPM --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFileNamePrefix /data2/${params.outdir}/${fastq_file.baseName}.

    bedtools bamtobed -i /data2/${params.outdir}/${fastq_file.baseName}.Aligned.toTranscriptome.out.bam > /data2/${params.outdir}/${fastq_file.baseName}.bed

    samtools sort /data2/${params.outdir}/${fastq_file.baseName}.Aligned.toTranscriptome.out.bam > /data2/${params.outdir}/${fastq_file.baseName}.Aligned.toTranscriptome.out.sorted.bam
    samtools index /data2/${params.outdir}/${fastq_file.baseName}.Aligned.toTranscriptome.out.sorted.bam

    """

}



process RunCTK {

  container 'wilfriedguiblet/ctk:v0.1'

  input:
    val fastq_file

  output:
    val "${fastq_file.baseName}"

  script:
    """
    cd /data2/
    export PERL5LIB=/opt/conda/lib/czplib
    /opt/conda/lib/ctk/tag2peak.pl \
      -big -ss \
      -p 0.05 --multi-test\
      --valley-seeking \
      --valley-depth 0.9 \
      /data2/${params.outdir}/${fastq_file.baseName}.bed /data2/${params.outdir}/${fastq_file.baseName}.peak.bed \
      --out-boundary /data2/${params.outdir}/${fastq_file.baseName}.uniq.peak.boundary.bed \
      --out-half-PH /data2/${params.outdir}/${fastq_file.baseName}.uniq.peak.halfPH.bed \
      --multi-test
    """
}


//process diffbind {

//  container 'wilfriedguiblet/dteg:v0.3'

//  input:
//    val fastq_file

//  output:
//    val "${fastq_file.baseName}"



workflow {

  RunRiboSeq(fastq_files)
  RunCTK(fastq_files)
}
