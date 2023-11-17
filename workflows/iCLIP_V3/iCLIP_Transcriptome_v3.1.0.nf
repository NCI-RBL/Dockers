"""
Author : Wilfried Guiblet. Blame him if it fails.
Update of : https://github.com/NCI-RBL/iCLIP

* Overview *
- Add-on to iCLIP pipeline
- Call peaks on transcriptome alignment
- Intersect with potential spliced-sites
- Build report

"""

// Necessary for syntax
nextflow.enable.dsl=2

// Create channels from the input paths

samplefiles_ch = Channel.fromList(params.samplenames)//.view { "value: $it" }
contrasts_ch = Channel.fromList(params.contrasts)//.view { "value: $it" }

// Create a channel with a unique value. Useful for processes that do not iterate through multiple samples.
unique_ch = Channel.fromList(['unique'])

// ************* The following parameters are imported from the file: nextflow.parameters.yaml *************

// determine which umi separator to use
if(params.multiplexflag == 'Y') {
    // demultiplexing ades rbc: to all demux files;
    params.umi_sep = "rbc:"}
else{
    // external demux uses an _
    params.umi_sep = params.umiSeparator}

params.manorm_w = params.MANormWidth
params.manorm_d = params.MNormDistance

// ************* End of parameter importation *************

process DeDup {
    """
    Sort and discard duplicates
    """
    container 'wilfriedguiblet/iclip:v3.1.0' // Use a Docker container

    input:
        val samplefile

    output:
        val samplefile

    shell:
        """
        set -exo pipefail

        # sort file
        samtools sort -m 80G -T !{params.workdir}/01_preprocess/02_alignment/ !{params.workdir}/01_preprocess/02_alignment/!{samplefile}_Aligned.toTranscriptome.out.bam -o !{params.workdir}/01_preprocess/02_alignment/!{samplefile}_Aligned.toTranscriptome.sortedByCoord.out.bam
        samtools index -@ !{params.featureCounts.threads} !{params.workdir}/01_preprocess/02_alignment/!{samplefile}_Aligned.toTranscriptome.sortedByCoord.out.bam

        # Run UMI Tools Deduplication
        echo "Using the following UMI seperator: !{params.umi_sep}"
        umi_tools dedup \\
            -I !{params.workdir}/01_preprocess/02_alignment/!{samplefile}_Aligned.toTranscriptome.sortedByCoord.out.bam \\
            --method unique \\
            --multimapping-detection-method=NH \\
            --umi-separator=!{params.umi_sep} \\
            -S !{params.workdir}/temp/!{samplefile}.toTranscriptome.unmasked.bam \\
            --log2stderr;
        
        # Sort and Index
        samtools sort --threads !{params.featureCounts.threads} -m 10G -T !{params.workdir}/temp/ \\
            !{params.workdir}/temp/!{samplefile}.toTranscriptome.unmasked.bam \\
            -o !{params.workdir}/02_bam/02_dedup/!{samplefile}.toTranscriptome.dedup.si.bam;
        samtools index -@ !{params.featureCounts.threads} !{params.workdir}/02_bam/02_dedup/!{samplefile}.toTranscriptome.dedup.si.bam;
        bedtools bamtobed -i !{params.workdir}/02_bam/02_dedup/!{samplefile}.toTranscriptome.dedup.si.bam > !{params.workdir}/03_peaks/01_bed/!{samplefile}.toTranscriptome.bed
        """
}


process CTK_Peak_Calling {
    """
    Peak calling using CTK.
    """

    container 'wilfriedguiblet/ctk:v0.1'

    input:
        val samplefile

    output:
        val samplefile

    shell:

        """
        export PERL5LIB=/opt/conda/lib/czplib

        /opt/conda/lib/ctk/tag2peak.pl \
        -big -ss \
        -p 0.001 --multi-test\
        --valley-seeking \
        --valley-depth 0.9 \
        !{params.workdir}/03_peaks/01_bed/!{samplefile}.toTranscriptome.bed !{params.workdir}/03_peaks/01_bed/!{samplefile}.toTranscriptome.peaks.bed \
        --out-boundary !{params.workdir}/03_peaks/01_bed/!{samplefile}.toTranscriptome.peaks.boundary.bed \
        --out-half-PH !{params.workdir}/03_peaks/01_bed/!{samplefile}.toTranscriptome.peaks.halfPH.bed \
        --multi-test \
        -minPH !{params.CTK.minimum_peak_height_Transcriptome}
        """
}

process Create_Safs {
    """
    Reformat BED into SAF.
    """

    container 'wilfriedguiblet/iclip:v3.1.0' // Use a Docker container

    input:
        val samplefile

    output:
        val samplefile

    shell:
        """
        set -exo pipefail
        awk '{{OFS="\\t"; if (\$3-\$2 >= 20) print \$1":"\$2"-"\$3"_"\$6,\$1,\$2,\$3,\$6}}' !{params.workdir}/03_peaks/01_bed/!{samplefile}.toTranscriptome.peaks.boundary.bed > !{params.workdir}/03_peaks/02_SAF/!{samplefile}.toTranscriptome.saf
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

    container 'wilfriedguiblet/iclip:v3.1.0' // Use a Docker container

    input:
        val samplefile

    output:
        val samplefile

    shell:
        """
        set -exo pipefail
        # Run for allreadpeaks
        featureCounts -F SAF \\
            -a !{params.workdir}/03_peaks/02_SAF/!{samplefile}.toTranscriptome.saf \\
            -O \\
            -J \\
            --fraction \\
            --minOverlap 1 \\
            -s 1 \\
            -T !{params.featureCounts.threads} \\
            -o /!{params.workdir}/03_peaks/03_counts/!{samplefile}_ALLreadpeaks_uniqueCounts.toTranscriptome.txt \\
            !{params.workdir}/02_bam/02_dedup/!{samplefile}.toTranscriptome.dedup.si.bam;
        featureCounts -F SAF \\
            -a !{params.workdir}/03_peaks/02_SAF/!{samplefile}.toTranscriptome.saf \\
            -M \\
            -O \\
            -J \\
            --fraction \\
            --minOverlap 1 \\
            -s 1 \\
            -T !{params.featureCounts.threads} \\
            -o !{params.workdir}/03_peaks/03_counts/!{samplefile}_ALLreadpeaks_FracMMCounts.toTranscriptome.txt \\
            !{params.workdir}/02_bam/02_dedup/!{samplefile}.toTranscriptome.dedup.si.bam;
        featureCounts -F SAF \\
            -a !{params.workdir}/03_peaks/02_SAF/!{samplefile}.toTranscriptome.saf \\
            -M \\
            -O \\
            --minOverlap 1 \\
            -s 1 \\
            -T !{params.featureCounts.threads} \\
            -o !{params.workdir}/03_peaks/03_counts/!{samplefile}_ALLreadpeaks_totalCounts.toTranscriptome.txt \\
            !{params.workdir}/02_bam/02_dedup/!{samplefile}.toTranscriptome.dedup.si.bam;
        """
}

process CombineCounts {
    """
    Combining the different type of counts done in FeatureCounts
    """

    container 'wilfriedguiblet/iclip:v3.1.0' // Use a Docker container

    input:
        val samplefile

    output:
        val samplefile

    shell:
        """
        # Usage: script input1 input2 input3 output
        python !{params.sourcedir}/05_countmerger.py \\
                    --uniqueCountsFile !{params.workdir}/03_peaks/03_counts/!{samplefile}_ALLreadpeaks_uniqueCounts.toTranscriptome.txt \\
                    --FracMMCountsFile !{params.workdir}/03_peaks/03_counts/!{samplefile}_ALLreadpeaks_FracMMCounts.toTranscriptome.txt \\
                    --totalCountsFile !{params.workdir}/03_peaks/03_counts/!{samplefile}_ALLreadpeaks_totalCounts.toTranscriptome.txt \\
                    --outName !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.toTranscriptome.txt
        """
}



process Peak_Annotation {
    """
    Annotate peaks Splice sites and Gene IDs
    """

    container 'wilfriedguiblet/iclip:v3.1.0' // Use a Docker container

    input:
        val samplefile

    output:
        val(samplefile)

    shell:
        """
        awk -v OFS='\t' '(NR>1) {print \$2, \$3, \$4, \$1, 0, \$5, \$6, \$7, \$8, \$9 }' \\
          !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.toTranscriptome.txt \\
            > !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.toTranscriptome.bed

        bedtools intersect -wa -wb -loj \\
          -a !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.toTranscriptome.bed \\
          -b !{params."${params.reference}".gencodedir}/GENCODE_VM23_knownGene_!{params.reference}.splice_sites.bed !{params."${params.reference}".gencodedir}/GENCODE_VM23_knownGene_!{params.reference}.5UTR.bed !{params."${params.reference}".gencodedir}/GENCODE_VM23_knownGene_!{params.reference}.3UTR.bed \\
          -names Splice_Sites 5UTR 3UTR \\
            > !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.!{params.reference}.intersect.toTranscriptome.bed

        python !{params.sourcedir}/MergeGeneIDs.py \\
           !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_AllRegions.!{params.reference}.intersect.toTranscriptome.bed \\
           !{params."${params.reference}".gencodepath} \\
           !{params.workdir}/04_annotation/02_peaks/!{samplefile}_!{params.peakid}readPeaks_annotation_complete.toTranscriptome.csv


        """

}


process MANORM_analysis {

    container 'wilfriedguiblet/iclip:v3.1.0' // Use a Docker container

    input:
       tuple val(sample), val(background), val(dummy)

    output:
        tuple val(sample), val(background)

    shell:
        """
        manorm \\
            --p1 "!{params.workdir}/03_peaks/01_bed/!{sample}.toTranscriptome.peaks.boundary.bed" \\
            --p2 "!{params.workdir}/03_peaks/01_bed/!{background}.toTranscriptome.peaks.boundary.bed" \\
            --r1 "!{params.workdir}/03_peaks/01_bed/!{sample}.toTranscriptome.bed" \\
            --r2 "!{params.workdir}/03_peaks/01_bed/!{background}.toTranscriptome.bed" \\
            --s1 0 \\
            --s2 0 \\
            -p 1 \\
            -m 0 \\
            -w !{params.manorm_w} \\
            --summit-dis !{params.manorm_d} \\
            --wa \\
            -o !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome \\
            --name1 !{sample} \\
            --name2 !{background}

        awk -v OFS='\t' '{print \$1,\$2,\$3,\$4}' !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/output_filters/!{sample}_vs_!{background}_M_above_0.0_biased_peaks.bed > !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome.manorm.xls
        """    
}


process MANORM_Annotation {
    """
    Annotate MANORM peaks Splice sites and Gene IDs
    """

    container 'wilfriedguiblet/iclip:v3.1.0' // Use a Docker container

    input:
        tuple val(sample), val(background)

    output:
        tuple val(sample), val(background)

    shell:
        """

        awk '{{OFS="\\t"; if (\$3-\$2 >= 20) print \$1":"\$2"-"\$3,\$1,\$2,\$3,"+"}}' !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/output_filters/!{sample}_vs_!{background}_M_above_0.0_biased_peaks.bed > !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.saf

        featureCounts -F SAF \\
            -a !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.saf \\
            -O \\
            -J \\
            --fraction \\
            --minOverlap 1 \\
            -s 1 \\
            -T !{params.featureCounts.threads} \\
            -o !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.uniqueCounts.!{sample}.txt \\
            !{params.workdir}/02_bam/02_dedup/!{sample}.toTranscriptome.dedup.si.bam;

        featureCounts -F SAF \\
            -a !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.saf \\
            -M \\
            -O \\
            -J \\
            --fraction \\
            --minOverlap 1 \\
            -s 1 \\
            -T !{params.featureCounts.threads} \\
            -o !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.FracMMCounts.!{sample}.txt \\
            !{params.workdir}/02_bam/02_dedup/!{sample}.toTranscriptome.dedup.si.bam;

        featureCounts -F SAF \\
            -a !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.saf \\
            -M \\
            -O \\
            --minOverlap 1 \\
            -s 1 \\
            -T !{params.featureCounts.threads} \\
            -o !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.totalCounts.!{sample}.txt \\
            !{params.workdir}/02_bam/02_dedup/!{sample}.toTranscriptome.dedup.si.bam;

        featureCounts -F SAF \\
            -a !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.saf \\
            -O \\
            -J \\
            --fraction \\
            --minOverlap 1 \\
            -s 1 \\
            -T !{params.featureCounts.threads} \\
            -o !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.uniqueCounts.!{background}.txt \\
            !{params.workdir}/02_bam/02_dedup/!{background}.toTranscriptome.dedup.si.bam;

        featureCounts -F SAF \\
            -a !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.saf \\
            -M \\
            -O \\
            -J \\
            --fraction \\
            --minOverlap 1 \\
            -s 1 \\
            -T !{params.featureCounts.threads} \\
            -o !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.FracMMCounts.!{background}.txt \\
            !{params.workdir}/02_bam/02_dedup/!{background}.toTranscriptome.dedup.si.bam;

        featureCounts -F SAF \\
            -a !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.saf \\
            -M \\
            -O \\
            --minOverlap 1 \\
            -s 1 \\
            -T !{params.featureCounts.threads} \\
            -o !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.totalCounts.!{background}.txt \\
            !{params.workdir}/02_bam/02_dedup/!{background}.toTranscriptome.dedup.si.bam;


        python !{params.sourcedir}/05_countmerger.py \\
                    --uniqueCountsFile !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.uniqueCounts.!{sample}.txt \\
                    --FracMMCountsFile !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.FracMMCounts.!{sample}.txt \\
                    --totalCountsFile !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.totalCounts.!{sample}.txt \\
                    --outName !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.AllCounts.!{sample}.txt

        python !{params.sourcedir}/05_countmerger.py \\
                    --uniqueCountsFile !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.uniqueCounts.!{background}.txt \\
                    --FracMMCountsFile !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.FracMMCounts.!{background}.txt \\
                    --totalCountsFile !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.totalCounts.!{background}.txt \\
                    --outName !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.AllCounts.!{background}.txt


        awk -v OFS='\t' '(NR>1) {print \$2, \$3, \$4, \$1, 0, \$5, \$6, \$7, \$8, \$9 }' \\
          !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.AllCounts.!{sample}.txt \\
            > !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.AllCounts.!{sample}.bed

        awk -v OFS='\t' '(NR>1) {print \$2, \$3, \$4, \$1, 0, \$5, \$6, \$7, \$8, \$9 }' \\
          !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.AllCounts.!{background}.txt \\
            > !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.AllCounts.!{background}.bed


        bedtools intersect -wa -wb -loj \\
          -a !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.AllCounts.!{sample}.bed \\
          -b !{params."${params.reference}".gencodedir}/GENCODE_VM23_knownGene_!{params.reference}.splice_sites.bed !{params."${params.reference}".gencodedir}/GENCODE_VM23_knownGene_!{params.reference}.5UTR.bed !{params."${params.reference}".gencodedir}/GENCODE_VM23_knownGene_!{params.reference}.3UTR.bed \\
          -names Splice_Sites 5UTR 3UTR \\
            > !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.intersect.!{sample}.bed

        python !{params.sourcedir}/MergeGeneIDsMANORM.py \\
           !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.intersect.!{sample}.bed \\
           !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.peaks.AllCounts.!{background}.bed \\
           !{params."${params.reference}".gencodepath} \\
           !{params.workdir}/05_demethod/02_analysis/!{sample}_vs_!{background}.toTranscriptome/!{sample}_vs_!{background}.annotation_complete.toTranscriptome.csv

        #
        """

}




workflow {

    DeDup(samplefiles_ch) | CTK_Peak_Calling | Create_Safs | Feature_Counts | CombineCounts | Peak_Annotation
    MANORM_analysis(contrasts_ch.combine(Peak_Annotation.out.collect().toList())) | MANORM_Annotation


}