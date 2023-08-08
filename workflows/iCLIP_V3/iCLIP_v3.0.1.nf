"""
Author : Wilfried Guiblet. Blame him if it fails.
Update of : https://github.com/NCI-RBL/iCLIP

* Overview *
- Multiplexed samples are split based on provided barcodes and named using provide manifests, maximum 10 samples
- Adaptors are stripped from samples
- Samples are unzipped and split into smaller fastq files to increase speed
- Samples are aligned using NovaAlign
- SAM and BAM files are created
- Samples are merged
* Requirements *
- Read specific input requirements, and execution information on the Wikipage
located at: TBD

"""

// Necessary for syntax
nextflow.enable.dsl=2

// Create channels from the input paths

rawfiles_ch = Channel.fromList(params.rawfilesnames)//.view { "value: $it" }
samplefiles_ch = Channel.fromList(params.samplenames)//.view { "value: $it" }
contrasts_ch = Channel.fromList(params.contrasts)//.view { "value: $it" }

// Create a channel with a unique value. Useful for processes that do not iterate through multiple samples.
unique_ch = Channel.fromList(['unique'])



// ************* The following parameters are imported from the file: nextflow.parameters.yaml *************


//params.threads = '4' // Threads to use for multithreading. Use carefully. Needs to be transfered to yaml.

// Convert rRNA selection from Y/N to TRUE/FALSE
if (params.include_rRNA=="Y") {params.rrna_flag = "TRUE"}
else {params.rrna_flag = "FALSE"}


params.a_config = "${params.workdir}/config/annotation_config.txt"
params.manorm_w = params.MANormWidth
params.manorm_d = params.MNormDistance

if( params.min_reads_mapped > 1) {
    println "Count_threshold must be a decimal value, representing a percentage."
}

// determine which umi separator to use
if(params.multiplexflag == 'Y') {
    // demultiplexing ades rbc: to all demux files;
    params.umi_sep = "rbc:"}
else{
    // external demux uses an _
    params.umi_sep = params.umiSeparator}

// ************* End of parameter importation *************


process Create_Project_Annotations {
    """
    Generate annotation table once per project.
    Generate BED files of annotations.
    """

    container 'wilfriedguiblet/iclip:v3.0.3' // Use a Docker container

    input:
        val unique

    output:
        val unique

    shell:
        """
        Rscript !{params.sourcedir}/04_annotation.R \\
          --ref_species !{params.reference} \\
          --refseq_rRNA !{params.rrna_flag} \\
          --alias_path !{params."${params.reference}".aliaspath} \\
          --gencode_path !{params."${params.reference}".gencodepath} \\
          --refseq_path !{params."${params.reference}".refseqpath} \\
          --canonical_path !{params."${params.reference}".can_path} \\
          --intron_path !{params."${params.reference}".intronpath} \\
          --rmsk_path !{params."${params.reference}".rmskpath} \\
          --custom_path !{params."${params.reference}".additionalannopath} \\
          --out_dir !{params.workdir}/04_annotation/01_project/ \\
          --reftable_path !{params.sourcedir}/annotation_config.txt 

        awk -v OFS='\t' '(NR>1) {print \$6, \$7, \$8, \$11, \$12, \$10, \$13}' !{params."${params.reference}".rmskpath} \\
             > !{params.workdir}/04_annotation/01_project/rmsk.!{params.reference}.bed
        awk -v OFS='\t' '(NR>1) {print \$1, \$4, \$5, \$2, \$6, \$7, \$3, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15, \$16,\\
           \$17, \$18, \$19, \$20, \$21, \$22, \$23, \$24, \$25}' !{params."${params.reference}".gencodepath} \\
             > !{params.workdir}/04_annotation/01_project/gencode.!{params.reference}.bed
        cp !{params."${params.reference}".intronpath} !{params.workdir}/04_annotation/01_project/KnownGene_introns.!{params.reference}.bed
        """
}

process Init_ReadCounts_Reportfile {

    container 'wilfriedguiblet/iclip:v3.0.3' // Use a Docker container

    input:
        val unique

    output:
        val unique

    shell:
        """
        # create output file
        if [[ -f !{params.workdir}/00_QC/02_SamStats/qc_read_count_raw_values.txt ]]; then rm !{params.workdir}/00_QC/02_SamStats/qc_read_count_raw_values.txt ; fi 
        touch !{params.workdir}/00_QC/02_SamStats/qc_read_count_raw_values.txt
        """
}


process QC_Barcode {
    """
    Barcodes will be reviewed to ensure uniformtiy amongst samples.
    - generate counts of barcodes and output to text file
    - run python script that determines barcode expected and generates mismatches based on input
    - output barplot with top barcode counts

    --mpid clip3 must be changed to be a variable etracted from relevant manifest
    """

    container 'wilfriedguiblet/iclip:v3.0.3' // Use a Docker container

    input:
        tuple val(unique), val(rawfile)

    output:
        val rawfile

    shell:    
        """
        set -exo pipefail

        gunzip -c !{params.rawdir}/!{rawfile}.fastq.gz \\
            | awk 'NR%4==2 {{print substr(\$0, !{params.qc_barcode.start_pos}, !{params.qc_barcode.barcode_length});}}' \\
            | LC_ALL=C sort --buffer-size=!{params.qc_barcode.memory} --parallel=!{params.qc_barcode.threads} --temporary-directory='!{params.tempdir}' -n \\
            | uniq -c > !{params.workdir}/00_QC/01_Barcodes/!{rawfile}_barcode_counts.txt;        

        Rscript !{params.sourcedir}/02_barcode_qc.R \\
            --sample_manifest !{params.manifests.samples} \\
            --multiplex_manifest !{params.manifests.multiplex} \\
            --barcode_input !{params.workdir}/00_QC/01_Barcodes/!{rawfile}_barcode_counts.txt \\
            --mismatch !{params.mismatch} \\
            --mpid !{params.qc_barcode.mpid} \\
            --output_dir !{params.workdir}/00_QC/01_Barcodes/ \\
            --qc_dir !{params.workdir}/00_QC/01_Barcodes/
        """

}

process Demultiplex {
    """
    https://github.com/ulelab/ultraplex

    NOTE: our SLURM system does not allow the use of --sbatchcompression which is recommended
    for increase in speed with --ultra. When the --sbatchcompression is used on our system, files 
    do not get compressed and will be transferred using a significant amount of disc space. 

    file_name                   multiplex
    SIM_iCLIP_S1_R1_001.fastq   SIM_iCLIP_S1
    multiplex       sample          group       barcode     adaptor
    SIM_iCLIP_S1    Ro_Clip         CLIP        NNNTGGCNN   AGATCGGAAGAGCGGTTCAG
    SIM_iCLIP_S1    Control_Clip    CNTRL       NNNCGGANN   AGATCGGAAGAGCGGTTCAG
    """

    container 'wilfriedguiblet/iclip:v3.0.3' // Use a Docker container

    input:
        val rawfile

    output:
        val rawfile

    shell:
        """
        set -exo pipefail
        
        # run ultraplex to remove adaptors, separate barcodes
        # output files to tmp scratch dir
        ultraplex \\
            --threads !{params.demultiplex.threads} \\
            --barcodes !{params.manifests.barcode} \\
            --directory !{params.workdir}/01_preprocess/01_fastq/ \\
            --inputfastq !{params.rawdir}/!{rawfile}.fastq.gz \\
            --final_min_length !{params.demultiplex.filterlength} \\
            --phredquality !{params.demultiplex.phredQuality} \\
            --fiveprimemismatches !{params.mismatch} \\
            --ultra 
        """
}

process Star {
    """
    STAR Alignment
    https://github.com/alexdobin/STAR/releases

    """
    container 'wilfriedguiblet/iclip:v3.0.3' // Use a Docker container

    input:
        tuple val(rawfile), val(samplefile)

    output:
        val samplefile

    shell:
        """
        set -exo pipefail

        # STAR cannot handle sorting large files - allow samtools to sort output files
        STAR \\
        --runThreadN !{params.STAR.threads} \\
        --runMode alignReads \\
        --genomeDir !{params."${params.reference}".stardir} \\
        --sjdbGTFfile !{params."${params.reference}".stargtf} \\
        --readFilesCommand zcat \\
        --readFilesIn !{params.workdir}/01_preprocess/01_fastq/ultraplex_demux_!{samplefile}.fastq.gz \\
        --outFileNamePrefix !{params.workdir}/01_preprocess/02_alignment/!{samplefile}_ \\
        --outReadsUnmapped Fastx \\
        --outSAMtype BAM Unsorted \\
        --alignEndsType !{params.STAR.alignEndsType} \\
        --alignIntronMax !{params.STAR.alignIntronMax} \\
        --alignSJDBoverhangMin !{params.STAR.alignSJDBoverhangMin} \\
        --alignSJoverhangMin !{params.STAR.alignSJoverhangMin} \\
        --alignTranscriptsPerReadNmax !{params.STAR.alignTranscriptsPerReadNmax} \\
        --alignWindowsPerReadNmax !{params.STAR.alignWindowsPerReadNmax} \\
        --limitBAMsortRAM !{params.STAR.bamlimit} \\
        --limitOutSJcollapsed !{params.STAR.limitOutSJcollapsed} \\
        --outFilterMatchNmin !{params.STAR.outFilterMatchNmin} \\
        --outFilterMatchNminOverLread !{params.STAR.outFilterMatchNminOverLread} \\
        --outFilterMismatchNmax !{params.STAR.outFilterMismatchNmax} \\
        --outFilterMismatchNoverReadLmax !{params.STAR.outFilterMismatchNoverReadLmax} \\
        --outFilterMultimapNmax !{params.STAR.outFilterMultimapNmax} \\
        --outFilterMultimapScoreRange !{params.STAR.outFilterMultimapScoreRange} \\
        --outFilterScoreMin !{params.STAR.outFilterScoreMin} \\
        --outFilterType !{params.STAR.outFilterType} \\
        --outSAMattributes !{params.STAR.outSAMattributes} \\
        --outSAMunmapped !{params.STAR.outSAMunmapped} \\
        --outSJfilterCountTotalMin !{params.STAR.outSJfilterCountTotalMin.replace(",", " ")} \\
        --outSJfilterOverhangMin !{params.STAR.outSJfilterOverhangMin.replace(",", " ")} \\
        --outSJfilterReads !{params.STAR.outSJfilterReads} \\
        --seedMultimapNmax !{params.STAR.seedMultimapNmax} \\
        --seedNoneLociPerWindow !{params.STAR.seedNoneLociPerWindow} \\
        --seedPerReadNmax !{params.STAR.seedPerReadNmax} \\
        --seedPerWindowNmax !{params.STAR.seedPerWindowNmax} \\
        --sjdbScore !{params.STAR.sjdbScore} \\
        --winAnchorMultimapNmax !{params.STAR.winAnchorMultimapNmax} \\
        --quantMode !{params.STAR.quantmod}

        # sort file
        samtools sort -m 80G -T !{params.workdir}/01_preprocess/02_alignment/ !{params.workdir}/01_preprocess/02_alignment/!{samplefile}_Aligned.out.bam -o !{params.workdir}/01_preprocess/02_alignment/!{samplefile}_Aligned.sortedByCoord.out.bam

        # move final log file to output
        mv !{params.workdir}/01_preprocess/02_alignment/!{samplefile}_Log.final.out !{params.workdir}/log/STAR/!{samplefile}.log
        
        # move mates to unmapped file
        touch !{params.workdir}/01_preprocess/02_alignment/!{samplefile}.unmapped.out
        for f in !{params.workdir}/01_preprocess/02_alignment/!{samplefile}_Unmapped.out.mate*; do cat \$f >> !{params.workdir}/01_preprocess/02_alignment/!{samplefile}.unmapped.out; done
        """

}

process Index_Stats{
    """
    sort, index files
    run samstats on files
    """

    container 'wilfriedguiblet/iclip:v3.0.3' // Use a Docker container

    input:
        val samplefile

    output:
        val samplefile
    
    shell:
        """
        set -exo pipefail
        
        # Index
        cp !{params.workdir}/01_preprocess/02_alignment/!{samplefile}_Aligned.sortedByCoord.out.bam !{params.workdir}/02_bam/01_merged/!{samplefile}.si.bam
        samtools index -@ !{params.featureCounts.threads} !{params.workdir}/02_bam/01_merged/!{samplefile}.si.bam;
        
        # Run samstats
        samtools stats --threads !{params.featureCounts.threads} !{params.workdir}/02_bam/01_merged/!{samplefile}.si.bam > !{params.workdir}/00_QC/02_SamStats/!{samplefile}_samstats.txt
        """

}



process Check_ReadCounts {
    """
    In a recent project the incorrect species was selected and nearly 80% of all reads in all samples (N=6) were not mapped. 
    Rather than continuing with this type of potential low-quality sample, the pipeline should stop.

    http://www.htslib.org/doc/samtools-stats.html
    """

    container 'wilfriedguiblet/iclip:v3.0.3' // Use a Docker container

    input:
        val samplefile

    output:
        val samplefile

    shell:
        """
        for f in !{params.workdir}/00_QC/02_SamStats/!{samplefile}_samstats.txt; do
            # check samstats file to determine number of reads and reads mapped
            raw_count=`cat \$f | grep "raw total sequences" | awk -F"\t" '{{print \$3}}'`
            mapped_count=`cat \$f | grep "reads mapped:" | awk -F"\t" '{{print \$3}}'`
            found_percentage=\$((\$mapped_count / \$raw_count))

            # check the count against the set count_threshold, if counts found are lower than expected, fail
            fail=0
            if [ 1 -eq "\$(echo "\${{found_percentage}} < !{params.min_reads_mapped}" | bc)" ]; then
                flag="sample failed"
                fail=\$((fail + 1))
            else
                flag="sample passed"
            fi
            
            # put data into output
            echo "\$f\t\$found_percentage\t\$flag" >> !{params.workdir}/00_QC/02_SamStats/qc_read_count_raw_values.txt
        done

        # create output file
if [ 1 -eq "\$(echo "\${{fail}} > 0" | bc)" ]; then
            echo "Check sample log !{params.workdir}/00_QC/02_SamStats/qc_read_count_raw_values.txt to review what sample(s) failed" > !{params.workdir}/00_QC/02_SamStats/qc_read_count_check_fail.txt
        else
            touch !{params.workdir}/00_QC/02_SamStats/qc_read_count_check_pass.txt
        fi
        """


}

process FastQC {

    container 'wilfriedguiblet/iclip:v3.0.3' // Use a Docker container

    input:
        tuple val(rawfile), val(samplefile)

    output:
        val samplefile

    shell:
        """
         set -exo pipefail

        # run FASTQC
        fastqc !{params.workdir}/01_preprocess/01_fastq/ultraplex_demux_!{samplefile}.fastq.gz \\
               -o !{params.workdir}/00_QC/03_MultiQC/
        """
}

process QC_Screen_Validator {
    """
    #fastq screen
    - this will align first to human, mouse, bacteria then will align to rRNA
    must run fastq_screen as two separate commands - multiqc will merge values of rRNA with human/mouse
    http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/_build/html/index.html
    - fastq validator
    Quality-control step to ensure the input FastQC files are not corrupted or
    incomplete prior to running the entire workflow.
    @Input:
        Raw FastQ file (scatter)
    @Output:
        Log file containing any warnings or errors on file
    """

    container 'wilfriedguiblet/iclip:v3.0.3' // Use a Docker container

    input:
        val samplefile

    output:
        val "QC_Done"

    shell:
        """
        set -exo pipefail

        # Gzip input files
        gunzip -c !{params.workdir}/01_preprocess/01_fastq/ultraplex_demux_!{samplefile}.fastq.gz \\
         > !{params.workdir}/temp/!{samplefile}.fastq;
        
        # Run FastQ Screen
        fastq_screen --conf !{params.sourcedir}/fqscreen_species_config.conf \\
            --outdir !{params.workdir}/00_QC/04_QC_ScreenSpecies \\
            --threads !{params.QC_Screen_Validator.threads} \\
            --subset 1000000 \\
            --aligner bowtie2 \\
            --force \\
            !{params.workdir}/temp/!{samplefile}.fastq ;
        fastq_screen --conf !{params.sourcedir}/fqscreen_rrna_config.conf \\
            --outdir !{params.workdir}/00_QC/05_QC_ScreenRRNA \\
            --threads !{params.QC_Screen_Validator.threads} \\
            --subset 1000000 \\
            --aligner bowtie2 \\
            --force \\
            !{params.workdir}/temp/!{samplefile}.fastq ;
        
        # Remove tmp gzipped file
        rm !{params.workdir}/temp/!{samplefile}.fastq
        
        # Run FastQ Validator
        !{params.fastq_val} \\
            --disableSeqIDCheck \\
            --noeof \\
            --printableErrors 100000000 \\
            --baseComposition \\
            --avgQual \\
            --file !{params.workdir}/01_preprocess/01_fastq/ultraplex_demux_!{samplefile}.fastq.gz \\
                > !{params.workdir}/00_QC/!{samplefile}.validated.fastq.log ;
        """
}


process MultiQC {
    """
    merges FastQC reports for pre/post trimmed fastq files into MultiQC report
    https://multiqc.info/docs/#running-multiqc
    """

    container 'wilfriedguiblet/iclip:v3.0.3' // Use a Docker container

    input:
        val check

    shell:
        """
        set -exo pipefail

        multiqc -f -v \\
            -c !{params.sourcedir}/multiqc_config.yaml \\
            -d -dd 1 \\
            !{params.workdir}/00_QC/03_MultiQC \\
            !{params.workdir}/00_QC/05_QC_ScreenRRNA \\
            !{params.workdir}/00_QC/04_QC_ScreenSpecies \\
            -o !{params.workdir}/00_QC/

        """
}


// rule qc_troubleshoot:



process DeDup {
    """
    deduplicate reads
    sort,index dedup.bam file
    get header of dedup file
    """

    container 'wilfriedguiblet/iclip:v3.0.3' // Use a Docker container

    input:
        val samplefile

    output:
        val samplefile

    shell:
        """
        set -exo pipefail
 
        # Run UMI Tools Deduplication
        echo "Using the following UMI seperator: !{params.umi_sep}"
        umi_tools dedup \\
            -I !{params.workdir}/02_bam/01_merged/!{samplefile}.si.bam \\
            --method unique \\
            --multimapping-detection-method=NH \\
            --umi-separator=!{params.umi_sep} \\
            -S !{params.workdir}/temp/!{samplefile}.unmasked.bam \\
            --log2stderr;
        
        # Sort and Index
        samtools sort --threads !{params.featureCounts.threads} -m 10G -T !{params.workdir}/temp/ \\
            !{params.workdir}/temp/!{samplefile}.unmasked.bam \\
            -o !{params.workdir}/02_bam/02_dedup/!{samplefile}.dedup.si.bam;
        samtools index -@ !{params.featureCounts.threads} !{params.workdir}/02_bam/02_dedup/!{samplefile}.dedup.si.bam;
        """
}

process Remove_Spliced_Reads {
    """
    Remove spliced reads from genome-wide alignment.
    Spliced reads create spliced peaks and will be dealt with by mapping against the transcriptome.
    """

    container 'wilfriedguiblet/iclip:v3.0.3' // Use a Docker container

    input:
        val samplefile

    output:
        val samplefile

    shell:
        """
        samtools view -h !{params.workdir}/02_bam/02_dedup/!{samplefile}.dedup.si.bam | awk -v OFS="\t" '\$0 ~ /^@/{print \$0;next;} \$6 !~ /N/' | samtools view -b > !{params.workdir}/02_bam/02_dedup/!{samplefile}.filtered.bam
        samtools index -@ !{params.featureCounts.threads} !{params.workdir}/02_bam/02_dedup/!{samplefile}.filtered.bam
        """
}

process CTK_Peak_Calling {
    """
    Alternative peak calling using CTK.
    """

    container 'wilfriedguiblet/ctk:v0.1'

    input:
        val samplefile

    output:
        val samplefile

    shell:

        """
        export PERL5LIB=/opt/conda/lib/czplib
        bedtools bamtobed -i !{params.workdir}/02_bam/02_dedup/!{samplefile}.filtered.bam > /data2/03_peaks/01_bed/!{samplefile}.bed

        /opt/conda/lib/ctk/tag2peak.pl \
        -big -ss \
        -p 0.001 --multi-test\
        --valley-seeking \
        --valley-depth 0.9 \
        !{params.workdir}/03_peaks/01_bed/!{samplefile}.bed !{params.workdir}/03_peaks/01_bed/!{samplefile}.peaks.bed \
        --out-boundary !{params.workdir}/03_peaks/01_bed/!{samplefile}.peaks.boundary.bed \
        --out-half-PH !{params.workdir}/03_peaks/01_bed/!{samplefile}.peaks.halfPH.bed \
        --multi-test \
        -minPH {!params.CTK.minimum_peak_height}
        """
}

process Create_Safs {
    """
    Reformat BED into SAF.
    """

    container 'wilfriedguiblet/iclip:v3.0.3' // Use a Docker container

    input:
        val samplefile

    output:
        val samplefile

    shell:
        """
        set -exo pipefail
        awk '{{OFS="\\t"; print \$1":"\$2"-"\$3"_"\$6,\$1,\$2,\$3,\$6}}' !{params.workdir}/03_peaks/01_bed/!{samplefile}.peaks.boundary.bed > !{params.workdir}/03_peaks/02_SAF/!{samplefile}.saf
        """
}

process Feature_Counts {
    """
    Unique reads (fractional counts correctly count splice reads for each peak.
    When peaks counts are combined for peaks connected by splicing in Rscript)
    Include Multimap reads - MM reads given fractional count based on # of mapping
    locations. All spliced reads also get fractional count. So Unique reads can get
    fractional count when spliced peaks combined in R script the summed counts give
    whole count for the unique alignement in combined peak.
    http://manpages.ubuntu.com/manpages/bionic/man1/featureCounts.1.html
    Output summary
    - Differences within any folder (allreadpeaks or uniquereadpeaks) should ONLY be the counts column -
    as this represent the number of peaks that were uniquely identified (uniqueCounts) or the number of peaks MM (allFracMMCounts)
    - Differences within folders (03_allreadpeaks, 03_uniquereadpeaks) will be the peaks identified, as the first takes
    all reads as input and the second takes only unique reads as input
    """

    container 'wilfriedguiblet/iclip:v3.0.3' // Use a Docker container

    input:
        val samplefile

    output:
        val samplefile

    shell:
        """
        set -exo pipefail
        # Run for allreadpeaks
        featureCounts -F SAF \\
            -a !{params.workdir}/03_peaks/02_SAF/!{samplefile}.saf \\
            -O \\
            -J \\
            --fraction \\
            --minOverlap 1 \\
            -s 1 \\
            -T !{params.featureCounts.threads} \\
            -o /!{params.workdir}/03_peaks/03_counts/!{samplefile}_ALLreadpeaks_uniqueCounts.txt \\
            !{params.workdir}/02_bam/02_dedup/!{samplefile}.filtered.bam;
        featureCounts -F SAF \\
            -a !{params.workdir}/03_peaks/02_SAF/!{samplefile}.saf \\
            -M \\
            -O \\
            -J \\
            --fraction \\
            --minOverlap 1 \\
            -s 1 \\
            -T !{params.featureCounts.threads} \\
            -o !{params.workdir}/03_peaks/03_counts/!{samplefile}_ALLreadpeaks_FracMMCounts.txt \\
            !{params.workdir}/02_bam/02_dedup/!{samplefile}.filtered.bam;
        featureCounts -F SAF \\
            -a !{params.workdir}/03_peaks/02_SAF/!{samplefile}.saf \\
            -M \\
            -O \\
            --minOverlap 1 \\
            -s 1 \\
            -T !{params.featureCounts.threads} \\
            -o !{params.workdir}/03_peaks/03_counts/!{samplefile}_ALLreadpeaks_totalCounts.txt \\
            !{params.workdir}/02_bam/02_dedup/!{samplefile}.filtered.bam;
        """
}

process CombineCounts {
    """
    Combining the different type of counts done in FeatureCounts
    """

    container 'wilfriedguiblet/iclip:v3.0.3' // Use a Docker container

    input:
        val samplefile

    output:
        val samplefile

    shell:
        """
        # Usage: script input1 input2 input3 output
        python !{params.sourcedir}/05_countmerger.py \\
                    --uniqueCountsFile !{params.workdir}/03_peaks/03_counts/!{samplefile}_ALLreadpeaks_uniqueCounts.txt \\
                    --FracMMCountsFile !{params.workdir}/03_peaks/03_counts/!{samplefile}_ALLreadpeaks_FracMMCounts.txt \\
                    --totalCountsFile !{params.workdir}/03_peaks/03_counts/!{samplefile}_ALLreadpeaks_totalCounts.txt \\
                    --outName !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.txt
        """
}

process Peak_Annotation {
    """
    Annotate peaks with GeneCode, Introns, and RepeatMasker
    """

    container 'wilfriedguiblet/iclip:v3.0.3' // Use a Docker container

    input:
        val samplefile

    output:
        val(samplefile)

    shell:
        """
        awk -v OFS='\t' '(NR>1) {print \$2, \$3, \$4, \$1, 0, \$5, \$6, \$7, \$8, \$9 }' \\
          !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.txt \\
            > !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.bed

        bedtools intersect -s -wao \\
          -a !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.bed \\
          -b !{params.workdir}/04_annotation/01_project/rmsk.!{params.reference}.bed \\
            > !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.rmsk.!{params.reference}.intersect.SameStrand.bed

        bedtools intersect -s -wao \\
          -a !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.bed \\
          -b !{params.workdir}/04_annotation/01_project/gencode.!{params.reference}.bed \\
          | awk 'BEGIN {FS = "\t"; OFS = "\t"} \$14 != "." {split(\$14, arr, "_"); \$18 = arr[3]} \$14 == "." { \$18 = "." } 1' \\
            > !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.gencode.!{params.reference}.intersect.SameStrand.bed

        bedtools intersect -s -wao \\
          -a !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.bed \\
          -b !{params.workdir}/04_annotation/01_project/KnownGene_introns.!{params.reference}.bed \\
            > !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.KnownGene_introns.!{params.reference}.intersect.SameStrand.bed


        bedtools intersect -S -wao \\
          -a !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.bed \\
          -b !{params.workdir}/04_annotation/01_project/rmsk.!{params.reference}.bed \\
            > !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.rmsk.!{params.reference}.intersect.OppoStrand.bed

        bedtools intersect -S -wao \\
          -a !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.bed \\
          -b !{params.workdir}/04_annotation/01_project/gencode.!{params.reference}.bed \\
            > !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.gencode.!{params.reference}.intersect.OppoStrand.bed

        bedtools intersect -S -wao \\
          -a !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.bed \\
          -b !{params.workdir}/04_annotation/01_project/KnownGene_introns.!{params.reference}.bed \\
          | awk 'BEGIN {FS = "\t"; OFS = "\t"} \$14 != "." {split(\$14, arr, "_"); \$18 = arr[3]} \$14 == "." { \$18 = "." } 1' \\
            > !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.KnownGene_introns.!{params.reference}.intersect.OppoStrand.bed


        python !{params.sourcedir}/AnnotationFormat.py \\
          --SameStrandRMSK !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.rmsk.!{params.reference}.intersect.SameStrand.bed \\
          --SameStrandGenCode !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.gencode.!{params.reference}.intersect.SameStrand.bed \\
          --SameStrandIntrons !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.KnownGene_introns.!{params.reference}.intersect.SameStrand.bed \\
          --OppoStrandRMSK !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.rmsk.!{params.reference}.intersect.OppoStrand.bed \\
          --OppoStrandGenCode !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.gencode.!{params.reference}.intersect.OppoStrand.bed \\
          --OppoStrandIntrons !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.KnownGene_introns.!{params.reference}.intersect.OppoStrand.bed \\
          --Output !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_annotation_complete.txt

        """
}


process Annotation_Report {
    """
    generates an HTML report for peak annotations
    """

    //container 'wilfriedguiblet/iclip:v3.0.3' // Use a Docker container

    input:
        val samplefile

    output:
        val "AnnotationDone"

    shell:
        """
        echo 'boom'
        #Rscript -e 'library(rmarkdown); \
        #rmarkdown::render("!{params.workdir}/06_annotation.Rmd", \
        #    output_file = "!{params.workdir}/04_annotation/!{samplefile}_!{params.peakid}readPeaks_final_report.html", \
        #    params= list(samplename = "!{samplefile}", \
        #        peak_in = "!{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_annotation_complete.txt", \
        #        output_table = "!{params.workdir}/04_annotation/!{samplefile}_annotation_!{params.peakid}readPeaks_final_table.txt", \
        #        readdepth = "!{params.mincount}", \
        #        PeakIdnt = "!{params.peakid}"))'
        """        

}

process MANORM_analysis {

    container 'wilfriedguiblet/iclip:v3.0.3' // Use a Docker container

    input:
       tuple val(sample), val(background), val(dummy)

    output:
        tuple val(sample), val(background)

    shell:
        """
        manorm \\
            --p1 "!{params.workdir}/03_peaks/01_bed/!{sample}.peaks.boundary.bed" \\
            --p2 "!{params.workdir}/03_peaks/01_bed/!{background}.peaks.boundary.bed" \\
            --r1 "!{params.workdir}/03_peaks/01_bed/!{sample}.bed" \\
            --r2 "!{params.workdir}/03_peaks/01_bed/!{background}.bed" \\
            --s1 0 \\
            --s2 0 \\
            -p 1 \\
            -m 0 \\
            -w !{params.manorm_w} \\
            --summit-dis !{params.manorm_d} \\
            --wa \\
            -o !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background} \\
            --name1 !{sample} \\
            --name2 !{background}


        awk -v OFS='\t' '{print \$1,\$2,\$3,\$4}' !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}/output_filters/!{sample}_M_above_0.0_biased_peaks.bed > !{params.workdir}/05_demethod/02_analysis/!{sample}.manorm.bed

        
        # rename MANORM final output file
        #mv {params.base}{wildcards.group_id}_all_MAvalues.xls {output.mavals}

        # rename individual file names for each sample
        #mv {params.base}{params.gid_1}_MAvalues.xls {params.base}{params.gid_1}_{params.peak_id}readPeaks_MAvalues.xls
        #mv {params.base}{params.gid_2}_MAvalues.xls {params.base}{params.gid_2}_{params.peak_id}readPeaks_MAvalues.xls

        # mv folders of figures, filters, tracks to new location
        # remove folders if they already exist
        #if [[ -d {params.base}output_figures_{params.peak_id}readPeaks ]]; then rm -r {params.base}output_figures_{params.peak_id}readPeaks; fi
        #if [[ -d {params.base}output_filters_{params.peak_id}readPeaks ]]; then rm -r {params.base}output_filters_{params.peak_id}readPeaks; fi
        #if [[ -d {params.base}output_tracks_{params.peak_id}readPeaks ]]; then rm -r {params.base}output_tracks_{params.peak_id}readPeaks; fi
        #mv {params.base}output_figures {params.base}output_figures_{params.peak_id}readPeaks
        #mv {params.base}output_filters {params.base}output_filters_{params.peak_id}readPeaks
        #mv {params.base}output_tracks {params.base}output_tracks_{params.peak_id}readPeaks
        """    

}

process Manorm_Report {

    container 'wilfriedguiblet/iclip:v3.0.3' // Use a Docker container

    input:
       tuple val(sample), val(background)

    output:
        tuple val(sample), val(background)

    shell:
        """
        set -exo pipefail
        
        ### for sample vs Bg compairson
        featureCounts -F SAF \\
            -a out_dir,'03_peaks','02_SAF','{gid_1}_' + peak_id + 'readPeaks.SAF' \\
            -O \\
            --fraction \\
            --minOverlap 1 \\
            -s 1 \\
            -T {threads} \\
            -o {output.bkUniqcountsmplPk} \\
            {params.bkbam}
        featureCounts -F SAF \\
            -a {params.smplSAF} \\
            -M \\
            -O \\
            --fraction \\
            --minOverlap 1 \\
            -s 1 \\
            -T {threads} \\
            -o {output.bkMMcountsmplPk} \\
            {params.bkbam}
        Rscript {params.script} \\
            --samplename {params.gid_1} \\
            --background {params.gid_2} \\
            --peak_anno_g1 {params.anno_dir}/{params.gid_1}_annotation_{params.peak_id}readPeaks_final_table.txt \\
            --peak_anno_g2 {params.anno_dir}/{params.gid_2}_annotation_{params.peak_id}readPeaks_final_table.txt \\
            --Smplpeak_bkgroundCount_MM {output.bkMMcountsmplPk} \\
            --Smplpeak_bkgroundCount_unique {output.bkUniqcountsmplPk} \\
            --pos_manorm {params.de_dir}/{wildcards.group_id}/{wildcards.group_id}_P/{params.gid_1}_{params.peak_id}readPeaks_MAvalues.xls \\
            --neg_manorm {params.de_dir}/{wildcards.group_id}/{wildcards.group_id}_N/{params.gid_1}_{params.peak_id}readPeaks_MAvalues.xls \\
            --output_file {output.post_proc}
        
        ### for Bg vs sample compairson
        featureCounts -F SAF \\
            -a {params.bkSAF} \\
            -O \\
            --fraction \\
            --minOverlap 1 \\
            -s 1 \\
            -T {threads} \\
            -o {output.smplUniqcountbkPk} \\
            {params.smplbam}
        featureCounts -F SAF \\
            -a {params.bkSAF} \\
            -M \\
            -O \\
            --fraction \\
            --minOverlap 1 \\
            -s 1 \\
            -T {threads} \\
            -o {output.smplMMcountbkPk} \\
            {params.smplbam}
        Rscript {params.script} \\
            --samplename {params.gid_2} \\
            --background {params.gid_1} \\
            --peak_anno_g1 {params.anno_dir}/{params.gid_2}_annotation_{params.peak_id}readPeaks_final_table.txt \\
            --peak_anno_g2 {params.anno_dir}/{params.gid_1}_annotation_{params.peak_id}readPeaks_final_table.txt \\
            --Smplpeak_bkgroundCount_MM {output.smplMMcountbkPk} \\
            --Smplpeak_bkgroundCount_unique {output.smplUniqcountbkPk} \\
            --pos_manorm {params.de_dir}/{wildcards.group_id}/{wildcards.group_id}_P/{params.gid_2}_{params.peak_id}readPeaks_MAvalues.xls \\
            --neg_manorm {params.de_dir}/{wildcards.group_id}/{wildcards.group_id}_N/{params.gid_2}_{params.peak_id}readPeaks_MAvalues.xls \\
            --output_file {output.post_procRev}
        """

}


workflow {

    Create_Project_Annotations(unique_ch) | Init_ReadCounts_Reportfile

    rawfiles_tuple = Init_ReadCounts_Reportfile.out.combine(rawfiles_ch)
    QC_Barcode(rawfiles_tuple) | Demultiplex

    samplefiles_tuple = Demultiplex.out.combine(samplefiles_ch)
    Star(samplefiles_tuple) | Index_Stats | Check_ReadCounts | DeDup | Remove_Spliced_Reads | CTK_Peak_Calling | Create_Safs | Feature_Counts | CombineCounts

    FastQC(samplefiles_tuple) | QC_Screen_Validator 
    MultiQC(QC_Screen_Validator.out.unique())

    Peak_Annotation(CombineCounts.out) | Annotation_Report

    MANORM_analysis(contrasts_ch.combine(Peak_Annotation.out.collect().toList()))


}
