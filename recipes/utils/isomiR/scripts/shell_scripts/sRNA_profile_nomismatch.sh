
echo "Start processing all _ready files in ./read_files/ folder..."

mkdir -p analysis_files
mkdir -p analysis_results

cd ready_files

ls -m1 *_ready > samples

cd ../analysis_files

printf "\ttotal\trRNA\ttRNA\tsnoRNA\tmiRNA\tmRNA\tothers_ref\tmycoplasma_H\tunmappable\thc_miRNA\n" > ../analysis_results/small_RNA_profile.txt

cat ../ready_files/samples | while read A
do
    echo "processing $A..."
    echo "mapping to rRNA..mode -n 0 -l 20..."
    bowtie -n 0 -l 20  --best --norc /data/RBL_NCI/Gu/miRNA_NGS/index/human/rRNA   --un "$A"_rm_rRNA.fastq  --al  "$A"_rRNA.fastq       ../ready_files/"$A" temp.txt
    echo "mapping to tRNA..mode -n 0 -l 20..."
    bowtie -n 0 -l 20 --best --norc /data/RBL_NCI/Gu/miRNA_NGS/index/human/tRNA   --un "$A"_rm_rRNA_tRNA.fastq  --al  "$A"_tRNA.fastq       "$A"_rm_rRNA.fastq temp.txt
    echo "mapping to snoRNA..mode -n 0 -l 20..."
    bowtie -n 0 -l 20 --best --norc /data/RBL_NCI/Gu/miRNA_NGS/index/human/snoRNA --un "$A"_rm_rRNA_tRNA_snoRNA.fastq  --al  "$A"_snoRNA.fastq     "$A"_rm_rRNA_tRNA.fastq temp.txt
    echo "mapping to miRNA_precusor.. n mode mismatch 0 seed length 20..."
    bowtie -n 0 -l 20 --best --norc --sam /data/RBL_NCI/Gu/miRNA_NGS/index/miRNA_Aug_2018/hsa_hairpin_all_2018  --un "$A"_rm_rRNA_tRNA_snoRNA_miRNA.fastq  --al  "$A"_miRNA.fastq  "$A"_rm_rRNA_tRNA_snoRNA.fastq    "$A"_miRNA.sam
    bowtie -n 0 -l 20 --best --norc --sam /data/RBL_NCI/Gu/miRNA_NGS/index/miRNA_Aug_2018/hsa_hairpin_hc_2018   --al  "$A"_hc_miRNA.fastq  "$A"_miRNA.fastq    "$A"_hc_miRNA.sam

     echo "mapping to mRNA.. mode -n 0 -l 20..."
    bowtie -n 0 -l 20 --best --norc /data/RBL_NCI/Gu/miRNA_NGS/index/human/refMRNA          --un "$A"_rm_rRNA_tRNA_snoRNA_miRNA_mRNA.fastq --al "$A"_mRNA.fastq   "$A"_rm_rRNA_tRNA_snoRNA_miRNA.fastq temp.txt
    echo "mapping to refSeq.. mode -n 0 -l 20..."
    bowtie -n 0 -l 20 --best --norc /data/RBL_NCI/Gu/miRNA_NGS/index/human/ref_transcripts          --un "$A"_others.fastq --al "$A"_others_ref.fastq   "$A"_rm_rRNA_tRNA_snoRNA_miRNA_mRNA.fastq temp.txt
    echo "processing $A..."
    echo "mapping to mycoplasma H..mode v1..."
    bowtie -v 1 --sam /data/RBL_NCI/Gu/miRNA_NGS/index/Mycoplasma_hyo --un "$A"_unmappable.fastq  --al "$A"_mycoplasmaH.fastq   ./"$A"_others.fastq  "$A"_mycoplasma.sam



    total=`cat ../ready_files/"$A" | wc -l`
    rRNA=`cat "$A"_rRNA.fastq | wc -l`
    tRNA=`cat "$A"_tRNA.fastq | wc -l`
    snoRNA=`cat "$A"_snoRNA.fastq | wc -l`
    miRNA=`cat "$A"_miRNA.fastq | wc -l`
    mRNA=`cat "$A"_mRNA.fastq | wc -l`
    others_ref=`cat "$A"_others_ref.fastq | wc -l`
    mycoplasma=`cat "$A"_mycoplasmaH.fastq | wc -l`
    unmappable=`cat "$A"_unmappable.fastq | wc -l`
    hc_miRNA=`cat "$A"_hc_miRNA.fastq | wc -l`


    let total=total/4
    let rRNA=rRNA/4
    let tRNA=tRNA/4
    let snoRNA=snoRNA/4
    let miRNA=miRNA/4
    let mRNA=mRNA/4
    let others_ref=others_ref/4
    let mycoplasma=mycoplasma/4
    let unmappable=unmappable/4
    let hc_miRNA=hc_miRNA/4

    rm temp.txt
    rm "$A"_rm_rRNA.fastq
    rm "$A"_rm_rRNA_tRNA.fastq
    rm "$A"_rm_rRNA_tRNA_snoRNA.fastq
    rm "$A"_rm_rRNA_tRNA_snoRNA_miRNA.fastq
    rm "$A"_rm_rRNA_tRNA_snoRNA_miRNA_mRNA.fastq
    rm  "$A"_others.fastq



printf "%s\t%d\t%d\t%d\t%d\t%d\t%d     \t%d\t%d\t%d\t%d\n" $A $total $rRNA $tRNA $snoRNA $miRNA $mRNA $others_ref $mycoplasma $unmappable $hc_miRNA>> ../analysis_results/small_RNA_profile.txt

done

cd ..




