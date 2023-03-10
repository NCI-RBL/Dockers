
echo " Start counting miRNAs from high confidence sam files in ./analysis_files/"

mkdir miRNA_counts
cd miRNA_counts

    cat ../ready_files/samples | while read A
    do
       ~/Desktop/NGS/scripts/count_miRNA_from_SAM_v2.exe ~/Desktop/NGS/index/miRNA_Aug_2018/hairpin_hc_2018.fa ../analysis_files/"$A"_hc_miRNA.sam ~/Desktop/NGS/index/miRNA_Aug_2018/miRNA_combine_list.txt  > "$A"_hc_miRNA_count.txt
    done
cd ..



