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
//bamfiles = Channel.fromPath(params.bamfiles)
fastqfiles = Channel.fromPath('${params.workdir}/01_preprocess/01_fastq/*.fastq.gz')
bamfiles = Channel.fromPath('bamfiles/*bam')
bedfiles = Channel.fromPath('03_peaks/01_bed/*bed')

// Create a channel with a unique value. Useful for processes that do not iterate through multiple samples.
create_unique = Channel.fromList(['unique'])



// ************* The following parameters are imported from the file: nextflow.parameters.yaml *************


params.workdir = '' // Working directory. Needs to be transfered to yaml.
params.threads = 1 // Threads to use for multithreading. Use carefully. Needs to be transfered to yaml.
params.tempdir = './' // Temp directory. 

// Convert rRNA selection and splice junction selection from Y/N to TRUE/FALSE
if (params.include_rRNA=="Y") {params.rrna_flag = "TRUE"}
else {params.rrna_flag = "FALSE"}
if (params.splicejunction=="Y") {params.sp_junc = "TRUE"}
else {params.sp_junc = "FALSE"}
//

params.a_path =  params."${params.reference}".aliaspath
params.g_path =  params."${params.reference}".gencodepath
params.rs_path =  params."${params.reference}".refseqpath
params.c_path =  params."${params.reference}".can_path
params.i_path =  params."${params.reference}".intronpath
params.r_path =  params."${params.reference}".rmskpath
params.custom_path =  params."${params.reference}".additionalannopath

params.s_index =  params."${params.reference}".stardir
params.s_gtf =  params."${params.reference}".stargtf
params.s_atype = params.alignEndsType
params.s_intron = params.alignIntronMax
params.s_sjdb = params.alignSJDBoverhangMin
params.s_asj = params.alignSJoverhangMin
params.s_transc = params.alignTranscriptsPerReadNmax
params.s_windows = params.alignWindowsPerReadNmax
params.star_bam_limit = "50297600554"
params.s_sjcol = params.limitOutSJcollapsed
params.s_match = params.outFilterMatchNmin
params.s_readmatch = params.outFilterMatchNminOverLread
params.s_mismatch = params.outFilterMismatchNmax
params.s_readmm = params.outFilterMismatchNoverReadLmax
params.s_fmm = params.outFilterMultimapNmax
params.s_mmscore = params.outFilterMultimapScoreRange
params.s_score = params.outFilterScoreMin
params.s_ftype = params.outFilterType
params.s_att = params.outSAMattributes
params.s_unmap = params.outSAMunmapped
params.s_sjmin = params.outSJfilterCountTotalMin.replace(",", " ")
params.s_overhang = params.outSJfilterOverhangMin.replace(",", " ")
params.s_sjreads = params.outSJfilterReads
params.s_smm = params.seedMultimapNmax
params.s_loci = params.seedNoneLociPerWindow
params.s_read = params.seedPerReadNmax
params.s_wind = params.seedPerWindowNmax
params.s_sj = params.sjdbScore
params.s_anchor = params.winAnchorMultimapNmax
params.s_quantmod = params.quantmod

params.a_config = "${params.workdir}/config/annotation_config.txt"

// ************* End of parameter importation *************


process Create_Project_Annotations {
    """
    Generate annotation table once per project.
    """

    //container 'wilfriedguiblet/iclip:v3.0' // Use a Docker container

    input:
        val file

    //output:
        //val file

    shell:
        """
        Rscript !{params.workdir}/workflow/scripts/04_annotation.R \\
            --ref_species !{params.reference} \\
            --refseq_rRNA !{params.rrna_flag} \\
            --alias_path !{params.a_path} \\
            --gencode_path !{params.g_path} \\
            --refseq_path !{params.rs_path} \\
            --canonical_path !{params.c_path} \\
            --intron_path !{params.i_path} \\
            --rmsk_path !{params.r_path} \\
            --custom_path !{params.custom_path} \\
            --out_dir !{params.workdir}/04_annotation/01_project/ \\
            --reftable_path !{params.a_config}
        """
}



process QC_Barcode {
    """
    Barcodes will be reviewed to ensure uniformtiy amongst samples.
    - generate counts of barcodes and output to text file
    - run python script that determines barcode expected and generates mismatches based on input
    - output barplot with top barcode counts
    """

    input:
        val file

    output:
        val "${file.SimpleName}"

    shell:    

        """
        set -exo pipefail

        gunzip -c !{file.SimpleName}.fastq.gz \\
            | awk 'NR%4==2 {{print substr(\$0, !{params.qc_barcode.start_pos}, !{params.qc_barcode.barcode_length});}}' \\
            | LC_ALL=C sort --buffer-size=!{params.qc_barcode.memory} --parallel=!{params.qc_barcode.threads} --temporary-directory='!{params.tempdir} -n \\
            | uniq -c > !{params.workdir}/00_QC/01_Barcodes/!{file.SimpleName}_barcode_counts.txt;        
            
        Rscript !{params.workdir}/workflow/scripts/02_barcode_qc.R \\
            --sample_manifest !{params.manifests.samples} \\
            --multiplex_manifest !{params.manifests.multiplex} \\
            --barcode_input !{params.workdir}/00_QC/01_Barcodes/!{file.SimpleName}_barcode_counts.txt \\
            --mismatch !{params.mismatch} \\
            --mpid !{file.SimpleName} \\
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

    input:
        val file

    output:
        val file

    shell:
        """
        set -exo pipefail
        
        # run ultraplex to remove adaptors, separate barcodes
        # output files to tmp scratch dir
        ultraplex \\
            --threads !{params.demultiplex.threads} \\
            --barcodes {input.barcodes} \\
            --directory \$tmp_dir \\
            --inputfastq {input.f1} \\
            --final_min_length {params.ml} \\
            --phredquality {params.pq} \\
            --fiveprimemismatches {params.mm} \\
            --ultra 

        # move files to final location after they are zipped
        mv \$tmp_dir/* {params.out_dir}
        """

}



process Star {
    """
    STAR Alignment
    https://github.com/alexdobin/STAR/releases

    """
    ///container 'wilfriedguiblet/iclip:v3.0'

    input:
        val file

    output:
        val "${file.SimpleName}"


    shell:
        """
        set -exo pipefail


        # STAR cannot handle sorting large files - allow samtools to sort output files
        STAR \
        --runMode alignReads \
        --genomeDir !{params.s_index} \
        --sjdbGTFfile !{params.s_gtf} \
        --readFilesCommand zcat \
        --readFilesIn !{params.workdir}/01_preprocess/01_fastq/!{file}.fastq.gz' \
        --outFileNamePrefix !{params.workdir}/01_preprocess/07_rscripts/ \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM Unsorted \
        --alignEndsType !{params.s_atype} \
        --alignIntronMax !{params.s_intron} \
        --alignSJDBoverhangMin !{params.s_sjdb} \
        --alignSJoverhangMin !{params.s_asj} \
        --alignTranscriptsPerReadNmax !{params.s_transc} \
        --alignWindowsPerReadNmax !{params.s_windows} \
        --limitBAMsortRAM !{params.star_bam_limit} \
        --limitOutSJcollapsed !{params.s_sjcol} \
        --outFilterMatchNmin !{params.s_match} \
        --outFilterMatchNminOverLread !{params.s_readmatch} \
        --outFilterMismatchNmax !{params.s_mismatch} \
        --outFilterMismatchNoverReadLmax !{params.s_readmm} \
        --outFilterMultimapNmax !{params.s_fmm} \
        --outFilterMultimapScoreRange !{params.s_mmscore} \
        --outFilterScoreMin !{params.s_score} \
        --outFilterType !{params.s_ftype} \
        --outSAMattributes !{params.s_att} \
        --outSAMunmapped !{params.s_unmap} \
        --outSJfilterCountTotalMin !{params.s_sjmin} \
        --outSJfilterOverhangMin !{params.s_overhang} \
        --outSJfilterReads !{params.s_sjreads} \
        --seedMultimapNmax !{params.s_smm} \
        --seedNoneLociPerWindow !{params.s_loci} \
        --seedPerReadNmax !{params.s_read} \
        --seedPerWindowNmax !{params.s_wind} \
        --sjdbScore !{params.s_sj} \
        --winAnchorMultimapNmax !{params.s_anchor} \
        --quantMode !{params.s_quantmod}

        # sort file
        samtools sort -m 80G -T !{params.workdir}/01_preprocess/07_rscripts/ !{params.workdir}/01_preprocess/07_rscripts/!{file}_Aligned.out.bam -o !{params.workdir}/01_preprocess/07_rscripts/!{file}_Aligned.sortedByCoord.out.bam

        # move STAR files and final log file to output
        mv !{params.workdir}/01_preprocess/07_rscripts/!{file}_Aligned.sortedByCoord.out.bam !{params.workdir}/01_preprocess/02_alignment/!{file}_Aligned.sortedByCoord.out.bam
        mv !{params.workdir}/01_preprocess/07_rscripts/!{file}_Log.final.out !{params.workdir}/log/STAR/!{file}.log'
        
        # move mates to unmapped file
        touch !{params.workdir}/01_preprocess/02_alignment/01_unmapped/!{file}.unmapped.out
        for f in !{params.workdir}/01_preprocess/07_rscripts/!{file}_Unmapped.out.mate*; do cat $f >> !{params.workdir}/01_preprocess/02_alignment/01_unmapped/!{file}.unmapped.out; done
        """

}


process CTK_Peak_Calling {
    """
    Alternative peak calling using CTK.
    """

    container 'wilfriedguiblet/iclip:v3.0'

    input:
        val file

    output:
        val file

    shell:

        """
        export PERL5LIB=/opt/conda/lib/czplib
        bedtools bamtobed -i /data2/bamfiles/!{file}.bam > /data2/bedfiles/!{file}.bed

        /opt/conda/lib/ctk/tag2peak.pl \
        -big -ss \
        -p 0.05 --multi-test\
        --valley-seeking \
        --valley-depth 0.9 \
        /data2/bedfiles/!{file}.bed /data2/peaks/!{file}.peaks.bed \
        --out-boundary /data2/peaks/!{file}.peaks.boundary.bed \
        --out-half-PH /data2/peaks/!{file}.peaks.halfPH.bed \
        --multi-test
        """

}



process Create_Safs {
    """
    Reformat BED into SAF.
    """

    input:
        val file

    output:
        val "${file.SimpleName}"

    shell:
        """
        set -exo pipefail
        awk '{{OFS="\\t"; print \$1":"\$2"-"\$3"_"\$6,\$1,\$2,\$3,\$6}}' !{file} > !{params.workdir}/03_peaks/02_SAF/!{file.SimpleName}.saf
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

    container 'wilfriedguiblet/iclip:v3.0'

    input:
        val file

    output:
        val file

    shell:
        """
        set -exo pipefail
        # Run for allreadpeaks
        featureCounts -F SAF \\
            -a /data2/03_peaks/02_SAF/!{file}.saf \\
            -O \\
            -J \\
            --fraction \\
            --minOverlap 1 \\
            -s 1 \\
            -T !{params.threads} \\
            -o /data2/03_peaks/03_counts/!{file}_ALLreadpeaks_uniqueCounts.txt \\
            /data2/bamfiles/!{file}.dedup.si.bam;
        featureCounts -F SAF \\
            -a /data2/03_peaks/02_SAF/!{file}.saf \\
            -M \\
            -O \\
            -J \\
            --fraction \\
            --minOverlap 1 \\
            -s 1 \\
            -T !{params.threads} \\
            -o /data2/03_peaks/03_counts/!{file}_ALLreadpeaks_FracMMCounts.txt \\
            /data2/bamfiles/!{file}.dedup.si.bam;
        featureCounts -F SAF \\
            -a /data2/03_peaks/02_SAF/!{file}.saf \\
            -M \\
            -O \\
            --minOverlap 1 \\
            -s 1 \\
            -T !{params.threads} \\
            -o /data2/03_peaks/03_counts/!{file}_ALLreadpeaks_totalCounts.txt \\
            /data2/bamfiles/!{file}.dedup.si.bam;
        # Run for uniquereadpeaks
        featureCounts -F SAF \\
            -a /data2/03_peaks/02_SAF/!{file}.saf \\
            -O \\
            -J \\
            --fraction \\
            --minOverlap 1 \\
            -s 1 \\
            -T !{params.threads} \\
            -o /data2/03_peaks/03_counts/!{file}_UNIQUEreadpeaks_uniqueCounts.txt \\
            /data2/bamfiles/!{file}.dedup.si.bam;
        featureCounts -F SAF \\
            -a /data2/03_peaks/02_SAF/!{file}.saf\\
            -M \\
            -O \\
            -J \\
            --fraction \\
            --minOverlap 1 \\
            -s 1 \\
            -T !{params.threads} \\
            -o /data2/03_peaks/03_counts/!{file}_UNIQUEreadpeaks_FracMMCounts.txt \\
            /data2/bamfiles/!{file}.dedup.si.bam;
        featureCounts -F SAF \\
            -a /data2/03_peaks/02_SAF/!{file}.saf \\
            -M \\
            -O \\
            --minOverlap 1 \\
            -s 1 \\
            -T !{params.threads} \\
            -o /data2/03_peaks/03_counts/!{file}_UNIQUEreadpeaks_totalCounts.txt \\
            /data2/bamfiles/!{file}.dedup.si.bam   
        """
}

process Peak_Junction {
    """
    find peak junctions, annotations peaks, merges junction and annotation information
    """

    input:
        val file

    output:
        val file

    shell:
        """
        #bash script to run bedtools and get site2peak lookuptable
        bash !{params.workdir}/workflow/scripts/05_get_site2peak_lookup.sh \\
             !{params.workdir}/03_peaks/03_counts/!{file}_!{params.peakid}readpeaks_FracMMCounts.txt.jcounts \\
             !{params.workdir}/03_peaks/03_counts/!{file}_!{params.peakid}readpeaks_FracMMCounts.txt \\
             !{file}_!{params.peakid} \\
             !{params.workdir}/04_annotation/02_peaks/ \\
             !{params.workdir}/workflow/scripts/05_jcounts2peakconnections.py

        # above bash script will create {output.splice_table}
        Rscript !{params.workdir}/workflow/scripts/05_Anno_junctions.R \\
            --rscript !{params.workdir}/workflow/scripts/05_peak_annotation_functions.R \\
            #--rscript !{params.workdir}/workflow/scripts/05_peak_annotation_functions_V2.3.R \\
            --peak_type !{params.peakid} \\
            --peak_unique !{params.workdir}/03_peaks/03_counts/!{file}_!{params.peakid}readpeaks_uniqueCounts.txt \\
            --peak_all !{params.workdir}/03_peaks/03_counts/!{file}_!{params.peakid}readpeaks_FracMMCounts.txt \\
            --peak_total !{params.workdir}/03_peaks/03_counts/!{file}_!{params.peakid}readpeaks_totalCounts.txt \\
            --join_junction !{params.sp_junc} \\
            --anno_anchor !{params.AnnoAnchor} \\
            --read_depth !{params.mincount} \\
            --demethod !{params.DEmethod} \\
            --sample_id !{file} \\
            --ref_species !{params.reference} \\
            --splice_table !{params.workdir}/04_annotation/02_peaks/!{file}_!{params.peakid}_connected_peaks.txt \\
            --tmp_dir !{params.workdir}/01_preprocess/07_rscripts/ \\
            --out_dir !{params.workdir}/04_annotation/02_peaks/ \\
            --out_file !{params.workdir}/04_annotation/02_peaks/!{file}_!{params.peakid}readPeaks_AllRegions.txt \\
            --out_dir_DEP !{params.workdir}/05_demethod/01_input/ \\
            --output_file_error !{params.workdir}/04_annotation/read_depth_error.txt
        """

}

process Peak_Transcripts {
    """
    find peak junctions, annotations peaks, merges junction and annotation information
    why is this the same description as Peak_Junction ?
    """

    input:
        val file

    output:
        val file

    shell:
        """
        Rscript !{params.workdir}/workflow/scripts/05_Anno_Transcript.R \\
            --rscript !{params.workdir}/workflow/scripts/05_peak_annotation_functions.R \\
            #--rscript !{params.workdir}/workflow/scripts/05_peak_annotation_functions_V2.3.R \\
            --peak_type !{params.peakid} \\
            --anno_anchor !{params.AnnoAnchor} \\
            --read_depth !{params.mincount} \\
            --sample_id !{file} \\
            --ref_species !{params.reference} \\
            --anno_dir !{params.workdir}/04_annotation/01_project/ \\
            --reftable_path !{params.a_config} \\
            --gencode_path !{params.g_path} \\
            --intron_path !{params.i_path} \\
            --rmsk_path !{params.r_path} \\
            --tmp_dir !{params.workdir}/01_preprocess/07_rscripts/ \\
            --out_dir !{params.workdir}/04_annotation/02_peaks/ \\
            --out_file !{params.workdir}/04_annotation/02_peaks/!{file}_!{params.peakid}readPeaks_AllRegions_transcripts_SameStrand.txt \\
            --anno_strand "SameStrand"

        Rscript !{params.workdir}/workflow/scripts/05_Anno_Transcript.R \\
            --rscript !{params.workdir}/workflow/scripts/05_peak_annotation_functions.R \\
            #--rscript !{params.workdir}/workflow/scripts/05_peak_annotation_functions_V2.3.R \\
            --peak_type !{params.peakid} \\
            --anno_anchor !{params.AnnoAnchor} \\
            --read_depth !{params.mincount} \\
            --sample_id !{file} \\
            --ref_species !{params.reference} \\
            --anno_dir !{params.workdir}/04_annotation/01_project/ \\
            --reftable_path !{params.a_config} \\
            --gencode_path !{params.g_path} \\
            --intron_path !{params.i_path} \\
            --rmsk_path !{params.r_path} \\
            --tmp_dir !{params.workdir}/01_preprocess/07_rscripts/ \\
            --out_dir !{params.workdir}/04_annotation/02_peaks/ \\
            --out_file !{params.workdir}/04_annotation/02_peaks/!{file}_!{params.peakid}readPeaks_AllRegions_transcripts_OppoStrand.txt \\
            --anno_strand "OppoStrand"
        """

}

process Peak_ExonIntron {
    """
    find peak junctions, annotations peaks, merges junction and annotation information
    why is this the same description as Peak_Junction ?
    """

    input:
        val file

    output:
        val file

    shell:
        """
        Rscript !{params.workdir}/workflow/scripts/05_Anno_ExonIntron.R \\
            --rscript !{params.workdir}/workflow/scripts/05_peak_annotation_functions.R \\
            --peak_type !{params.peakid} \\
            --anno_anchor !{params.AnnoAnchor} \\
            --read_depth !{params.mincount} \\
            --sample_id !{file} \\
            --ref_species !{params.reference} \\
            --anno_dir !{params.workdir}/04_annotation/01_project/ \\
            --reftable_path !{params.a_config} \\
            --gencode_path !{params.g_path} \\
            --intron_path !{params.i_path} \\
            --rmsk_path !{params.r_path} \\
            --tmp_dir !{params.workdir}/01_preprocess/07_rscripts/same \\
            --out_dir !{params.workdir}/04_annotation/02_peaks/ \\
            --out_file !{params.workdir}/04_annotation/02_peaks/!{file}_!{params.peakid}readPeaks_AllRegions_IntronExon_SameStrand.txt \\
            --anno_strand "SameStrand" 

        Rscript !{params.workdir}/workflow/scripts/05_Anno_ExonIntron.R \\
            --rscript !{params.workdir}/workflow/scripts/05_peak_annotation_functions.R \\
            --peak_type !{params.peakid} \\
            --anno_anchor !{params.AnnoAnchor} \\
            --read_depth !{params.mincount} \\
            --sample_id !{file} \\
            --ref_species !{params.reference} \\
            --anno_dir !{params.workdir}/04_annotation/01_project/ \\
            --reftable_path !{params.a_config} \\
            --gencode_path !{params.g_path} \\
            --intron_path !{params.i_path} \\
            --rmsk_path !{params.r_path} \\
            --tmp_dir !{params.workdir}/01_preprocess/07_rscripts/oppo \\
            --out_dir !{params.workdir}/04_annotation/02_peaks/ \\
            --out_file !{params.workdir}/04_annotation/02_peaks/!{file}_!{params.peakid}readPeaks_AllRegions_IntronExon_OppoStrand.txt \\
            --anno_strand "OppoStrand"
        """
}


process Peak_RMSK {
    """
    find peak junctions, annotations peaks, merges junction and annotation information
    why is this the same description as Peak_Junction ?
    """

    input:
        val file

    output:
        val file

    shell:
        """
        Rscript !{params.workdir}/workflow/scripts/05_Anno_RMSK.R \\
            --rscript !{params.workdir}/workflow/scripts/05_peak_annotation_functions.R \\
            --peak_type !{params.peakid} \\
            --anno_anchor !{params.AnnoAnchor} \\
            --read_depth !{params.mincount} \\
            --sample_id !{file} \\
            --ref_species !{params.reference} \\
            --anno_dir !{params.workdir}/04_annotation/01_project/ \\
            --reftable_path !{params.a_config} \\
            --gencode_path !{params.g_path} \\
            --intron_path !{params.i_path} \\
            --rmsk_path !{params.r_path} \\
            --tmp_dir !{params.workdir}/01_preprocess/07_rscripts/ \\
            --out_dir !{params.workdir}/04_annotation/02_peaks/ \\
            --out_file !{params.workdir}/04_annotation/02_peaks/!{file}_!{params.peakid}readPeaks_AllRegions_RMSK_SameStrand.txt \\
            --anno_strand "SameStrand"

        Rscript !{params.workdir}/workflow/scripts/05_Anno_RMSK.R \\
            --rscript !{params.workdir}/workflow/scripts/05_peak_annotation_functions.R \\
            --peak_type !{params.peakid} \\
            --anno_anchor !{params.AnnoAnchor} \\
            --read_depth !{params.mincount} \\
            --sample_id !{file} \\
            --ref_species !{params.reference} \\
            --anno_dir !{params.workdir}/04_annotation/01_project/ \\
            --reftable_path !{params.a_config} \\
            --gencode_path !{params.g_path} \\
            --intron_path !{params.i_path} \\
            --rmsk_path !{params.r_path} \\
            --tmp_dir !{params.workdir}/01_preprocess/07_rscripts/ \\
            --out_dir !{params.workdir}/04_annotation/02_peaks/ \\
            --out_file !{params.workdir}/04_annotation/02_peaks/!{file}_!{params.peakid}readPeaks_AllRegions_RMSK_OppoStrand.txt \\
            --anno_strand "OppoStrand"
        """
}

process Peak_Process {
    """
    find peak junctions, annotations peaks, merges junction and annotation information
    why is this the same description as Peak_Junction ?
    """

    input:
        val file

    output:
        val file

    shell:
        """
        Rscript !{params.workdir}/workflow/scripts/05_Anno_Process.R \\
            --rscript !{params.workdir}/workflow/scripts/05_peak_annotation_functions.R \\
            --peak_type !{params.peakid} \\
            --anno_anchor !{params.AnnoAnchor} \\
            --read_depth !{params.mincount} \\
            --sample_id !{file} \\
            --ref_species !{params.reference} \\
            --tmp_dir !{params.workdir}/01_preprocess/07_rscripts/ \\
            --out_dir !{params.workdir}/04_annotation/02_peaks/ \\
            --out_file !{params.workdir}/04_annotation/02_peaks/!{file}_!{params.peakid}readPeaks_annotation_complete.txt
        """

}


process Annotation_Report {
    """
    generates an HTML report for peak annotations
    """

    input:
        val file

    output:
        val file

    shell:
        """
        Rscript -e 'library(rmarkdown); \
        rmarkdown::render("!{params.workdir}/workflow/scripts/06_annotation.Rmd",
            output_file = "!{params.workdir}/04_annotation/!{file.SimpleName}_!{params.peakid}readPeaks_final_report.html", \
            params= list(samplename = "!{file}", \
                peak_in = "!{params.workdir}/04_annotation/02_peaks/!{file.SimpleName}_!{params.peakid}readPeaks_annotation_complete.txt", \
                output_table = "!{params.workdir}/04_annotation/!{file.SimpleName}_annotation_!{params.peakid}readPeaks_final_table.txt", \
                readdepth = "!{params.mincount}", \
                PeakIdnt = "!{params.peakid}"))'
        """        


}


workflow {
    //Create_Project_Annotations(create_unique)
    //Star(bamfiles)
    //Create_Safs(bedfiles)
    //Feature_Counts(Create_Safs.out)
    //Peak_Junction(Feature_Counts.out)
    //Peak_Transcripts(Peak_Junction.out)
    //Peak_ExonIntron(Peak_Transcripts.out)
    //Peak_RMSK(Peak_ExonIntron.out)
    //Peak_Process(Peak_RMSK.out)
    //Annotation_Report(Peak_Process.out)
    Annotation_Report(bedfiles)
}
