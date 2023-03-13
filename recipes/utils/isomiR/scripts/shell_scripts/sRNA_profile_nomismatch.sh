
echo "Start processing all _ready files in ./ready_files/ folder..."

mkdir -p /data2/analysis_files
mkdir -p /data2/analysis_results

cd /data2/ready_files
ls -m1 *_ready > samples
cd ../analysis_files

printf "\ttotal\trRNA\ttRNA\tsnoRNA\tmiRNA\tmRNA\tothers_ref\tmycoplasma_H\tunmappable\thc_miRNA\n" > /data2/analysis_results/small_RNA_profile.txt

cat /data2/ready_files/samples | while read A

do
    echo "processing $A..."
    echo "mapping to rRNA..mode -n 0 -l 20..."
    bowtie -n 0 -l 20 --best --norc /data2/NGS/index/human/rRNA --un /data2/analysis_files/"$A"_rm_rRNA.fastq --al /data2/analysis_files/"$A"_rRNA.fastq /data2/ready_files/"$A" /data2/analysis_files/temp.txt
    echo "mapping to tRNA..mode -n 0 -l 20..."
    bowtie -n 0 -l 20 --best --norc /data2/NGS/index/human/tRNA --un /data2/analysis_files/"$A"_rm_rRNA_tRNA.fastq --al /data2/analysis_files/"$A"_tRNA.fastq /data2/analysis_files/"$A"_rm_rRNA.fastq /data2/analysis_files/temp.txt
    echo "mapping to snoRNA..mode -n 0 -l 20..."
    bowtie -n 0 -l 20 --best --norc /data2/NGS/index/human/snoRNA --un /data2/analysis_files/"$A"_rm_rRNA_tRNA_snoRNA.fastq --al /data2/analysis_files/"$A"_snoRNA.fastq /data2/analysis_files/"$A"_rm_rRNA_tRNA.fastq /data2/analysis_files/temp.txt
    echo "mapping to miRNA_precusor.. n mode mismatch 0 seed length 20..."
    bowtie -n 0 -l 20 --best --norc --sam /data2/NGS/index/miRNA_Aug_2018/hsa_hairpin_all_2018 --un /data2/analysis_files/"$A"_rm_rRNA_tRNA_snoRNA_miRNA.fastq --al /data2/analysis_files/"$A"_miRNA.fastq /data2/analysis_files/"$A"_rm_rRNA_tRNA_snoRNA.fastq /data2/analysis_files/"$A"_miRNA.sam
    bowtie -n 0 -l 20 --best --norc --sam /data2/NGS/index/miRNA_Aug_2018/hsa_hairpin_hc_2018 --al /data2/analysis_files/"$A"_hc_miRNA.fastq /data2/analysis_files/"$A"_miRNA.fastq /data2/analysis_files/"$A"_hc_miRNA.sam

     echo "mapping to mRNA.. mode -n 0 -l 20..."
    bowtie -n 0 -l 20 --best --norc /data2/NGS/index/human/refMrna --un /data2/analysis_files/"$A"_rm_rRNA_tRNA_snoRNA_miRNA_mRNA.fastq --al /data2/analysis_files/"$A"_mRNA.fastq /data2/analysis_files/"$A"_rm_rRNA_tRNA_snoRNA_miRNA.fastq /data2/analysis_files/temp.txt
    echo "mapping to refSeq.. mode -n 0 -l 20..."
    bowtie -n 0 -l 20 --best --norc /data2/NGS/index/human/ref_transcripts --un /data2/analysis_files/"$A"_others.fastq --al /data2/analysis_files/"$A"_others_ref.fastq /data2/analysis_files/"$A"_rm_rRNA_tRNA_snoRNA_miRNA_mRNA.fastq /data2/analysis_files/temp.txt
    echo "processing $A..."
    echo "mapping to mycoplasma H..mode v1..."
    bowtie -v 1 --sam /data2/NGS/index/Mycoplasma_hyo --un /data2/analysis_files/"$A"_unmappable.fastq  --al /data2/analysis_files/"$A"_mycoplasmaH.fastq /data2/analysis_files/"$A"_others.fastq  /data2/analysis_files/"$A"_mycoplasma.sam

    total=`cat /data2/ready_files/"$A" | wc -l`
    rRNA=`cat /data2/analysis_files/"$A"_rRNA.fastq | wc -l`
    tRNA=`cat /data2/analysis_files/"$A"_tRNA.fastq | wc -l`
    snoRNA=`cat /data2/analysis_files/"$A"_snoRNA.fastq | wc -l`
    miRNA=`cat /data2/analysis_files/"$A"_miRNA.fastq | wc -l`
    mRNA=`cat /data2/analysis_files/"$A"_mRNA.fastq | wc -l`
    others_ref=`cat /data2/analysis_files/"$A"_others_ref.fastq | wc -l`
    mycoplasma=`cat /data2/analysis_files/"$A"_mycoplasmaH.fastq | wc -l`
    unmappable=`cat /data2/analysis_files/"$A"_unmappable.fastq | wc -l`
    hc_miRNA=`cat /data2/analysis_files/"$A"_hc_miRNA.fastq | wc -l`

    let total=total/4
    let rRNA=rRNA/4
    let tRNA=tRNA/4
    let snoRNA=snoRNA/4
    let miRNA=miRNA/4
    let mRNA=mRNA/4
    let others_ref=others_ref/4
    let mycoplasma=mycoplasma/4
    let unmappable=unmappable/4

printf "%s\t%d\t%d\t%d\t%d\t%d\t%d     \t%d\t%d\t%d\t%d\n" $A $total $rRNA $tRNA $snoRNA $miRNA $mRNA $others_ref $mycoplasma $unmappable $hc_miRNA>> /data2/analysis_results/small_RNA_profile.txt

done

