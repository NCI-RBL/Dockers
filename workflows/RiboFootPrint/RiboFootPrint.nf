// Author : Wilfried Guiblet

nextflow.enable.dsl=2


fastq_files = Channel.fromList(params.samples)

subparts = Channel.fromList(['5utr', 'cds', '3utr'])

// Processes

process RunRiboSeq {

  container 'wilfriedguiblet/dteg:v0.3'

  input:
    val fastq_file

  output:
    val fastq_file

  shell:
    """
    #seqtk trimfq -b 4 /data2/${params.fastqdir}/${fastq_file}.fastq.gz  > /data2/${params.outdir}/${fastq_file}.trimmed
    #skewer -t 26 -x ${params.linker} -l 15 /data2/${params.outdir}/${fastq_file}.trimmed

    #STAR --runThreadN 26 --readFilesIn /data2/${params.outdir}/${fastq_file}-trimmed.fastq --genomeDir /data2/index/star_hg19_ncrna/ --outReadsUnmapped Fastx --outSAMtype None --outFilterMismatchNoverLmax 0.3 --outFileNamePrefix /data2/${params.outdir}/${fastq_file}.ncrnaout.

    #STAR --runThreadN 26 --readFilesIn /data2/${params.outdir}/${fastq_file}.ncrnaout.Unmapped.out.mate1 --limitBAMsortRAM 600000000000 --genomeDir /data2/index/hg19/ --outFilterType BySJout --outFilterIntronMotifs RemoveNoncanonicalUnannotated --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.1 --outWigType wiggle read1_5p --outWigNorm RPM --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFileNamePrefix /data2/${params.outdir}/${fastq_file}.

    #samtools sort /data2/${params.outdir}/${fastq_file}.Aligned.toTranscriptome.out.bam > /data2/${params.outdir}/${fastq_file}.hg19.transcriptome.bam 
    #samtools index /data2/${params.outdir}/${fastq_file}.hg19.transcriptome.bam
    #samtools view -c /data2/${params.outdir}/${fastq_file}.hg19.transcriptome.bam > /data2/${params.outdir}/${fastq_file}.readcount
    bedtools bamtobed -i /data2/${params.outdir}/${fastq_file}.hg19.transcriptome.bam > /data2/${params.outdir}/${fastq_file}.hg19.transcriptome.bed

    awk -v OFS='\t' '{if (\$3-\$2 <= !{params.MaxReadLength}) print \$0}' /data2/${params.outdir}/${fastq_file}.hg19.transcriptome.bed > /data2/${params.outdir}/${fastq_file}.hg19.transcriptome.filtered.bed
    """

}

process ReadStarts {

  input:
    val fastq_file

  output:
    val fastq_file

  shell:
    """
    python ${params.workdir}/ReadStarts_V2.py --Infile ${params.workdir}/${params.outdir}/${fastq_file}.hg19.transcriptome.filtered.bed --Outfile ${params.workdir}/${params.outdir}/${fastq_file}.hg19.transcriptome.starts --TranscriptomeFile ${params.workdir}/hg19.transcriptome.bed --Shift ${params.Shift}
    """
}

process Relative_Aggregate {

  input:
    val fastq_file

  output:
    val fastq_file

  shell:
    """
    python ${params.workdir}/Relative_Aggregate_V2.py --Infile ${params.workdir}/${params.outdir}/${fastq_file}.hg19.transcriptome.starts --Outfile ${params.workdir}/${params.outdir}/${fastq_file}.hg19.transcriptome.RelativeAggregate --SubPartsFile ${params.workdir}/subparts.5utr20bp.txt
    """
}

process Plot_Aggregates {

  input:
    val fastq_file

  output:
    val fastq_file

  shell:
    """
    Rscript ${params.workdir}/PlotAggregates_V2.r ${params.workdir}/${params.outdir}/${fastq_file}.hg19.transcriptome.RelativeAggregate ${params.workdir}/${params.outdir}/ ${fastq_file} ${params.workdir}/${params.outdir}/${fastq_file}.readcount
    """
}


workflow {

  //RunRiboSeq(fastq_files) | ReadStarts | Relative_Aggregate | Plot_Aggregates
  Plot_Aggregates(fastq_files)
}
