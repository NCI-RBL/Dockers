
// Author : Wilfried Guiblet
// Update of : https://github.com/NCI-RBL/iCLIP


nextflow.enable.dsl=2

// Create a channel from the input path
//bamfiles = Channel.fromPath(params.bamfiles)
bamfiles = Channel.fromPath('bamfiles/*bam')
bedfiles = Channel.fromPath('03_peaks/01_bed/*bed')

create_unique = Channel.fromList(['unique'])

params.workdir = ''
params.threads = 1


// convert rRNA selection
if (params.include_rRNA=="Y") {params.rrna_flag = "TRUE"}
else {params.rrna_flag = "FALSE"}
//

// convert splice junction selection
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
params.a_config = "${params.workdir}/config/annotation_config.txt"


process Create_Project_Annotations {
    """
    generate annotation table once per project
    """

    //container 'wilfriedguiblet/iclip:v3.0'

    input:
        val file

    output:
        val file

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


process Create_Safs {
    """
    create SAF files
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

workflow {
    //Create_Project_Annotations(create_unique)
    Create_Safs(bedfiles)
    Feature_Counts(Create_Safs.out)
    Peak_Junction(Feature_Counts.out)
    Peak_Transcripts(Peak_Junction.out)
    Peak_ExonIntron(Peak_Transcripts.out)
    Peak_RMSK(Peak_ExonIntron.out)
    Peak_Process(Peak_RMSK.out)
}

