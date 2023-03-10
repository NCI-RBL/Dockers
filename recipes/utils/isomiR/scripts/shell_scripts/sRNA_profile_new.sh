
echo "Start processing all _ready files in ./read_files/ folder..."

mkdir analysis_files
mkdir analysis_results

cd ready_files

ls -m1 *_ready > samples

cd ../analysis_files

printf "\ttotal\trRNA\ttRNA\tsnoRNA\tmiRNA\tmRNA\tothers\tunmappable\n" > ../analysis_results/small_RNA_profile.txt

cat ../ready_files/samples | while read A
do
    echo "processing $A..."
    echo "mapping to rRNA..mode -n 1 -l 20..."
    ~/Desktop/NGS/tools/bowtie-1.1.1/bowtie -n 1 -l 20  --best --norc ~/Desktop/NGS/index/human/rRNA   --un "$A"_rm_rRNA.fastq  --al  "$A"_rRNA.fastq       ../ready_files/"$A" temp.txt
    echo "mapping to tRNA..mode -n l -l 20..."
    ~/Desktop/NGS/tools/bowtie-1.1.1/bowtie -n 1 -l 20 --best --norc ~/Desktop/NGS/index/human/tRNA   --un "$A"_rm_rRNA_tRNA.fastq  --al  "$A"_tRNA.fastq       "$A"_rm_rRNA.fastq temp.txt
    echo "mapping to snoRNA..mode -n 1 -l 20..."
    ~/Desktop/NGS/tools/bowtie-1.1.1/bowtie -n 1 -l 20 --best --norc ~/Desktop/NGS/index/human/snoRNA --un "$A"_rm_rRNA_tRNA_snoRNA.fastq  --al  "$A"_snoRNA.fastq     "$A"_rm_rRNA_tRNA.fastq temp.txt
    echo "mapping to miRNA_precusor.. n mode mismatch 1 seed length 20..."
    ~/Desktop/NGS/tools/bowtie-1.1.1/bowtie -n 1 -l 20 --best --norc --sam ~/Desktop/NGS/index/hsa_hairpin_all  --un "$A"_rm_rRNA_tRNA_snoRNA_miRNA.fastq  --al  "$A"_miRNA.fastq  "$A"_rm_rRNA_tRNA_snoRNA.fastq    "$A"_miRNA.sam
     echo "mapping to mRNA.. mode -n 1 -l 20..."
    ~/Desktop/NGS/tools/bowtie-1.1.1/bowtie -n 1 -l 20 --best --norc ~/Desktop/NGS/index/human/refMRNA          --un "$A"_rm_rRNA_tRNA_snoRNA_miRNA_mRNA.fastq --al "$A"_mRNA.fastq   "$A"_rm_rRNA_tRNA_snoRNA_miRNA.fastq temp.txt
    echo "mapping to refSeq.. mode -n 1 -l 20..."
    ~/Desktop/NGS/tools/bowtie-1.1.1/bowtie -n 1 -l 20 --best --norc ~/Desktop/NGS/index/human/ref_transcripts          --un "$A"_unmappable.fastq --al "$A"_others.fastq   "$A"_rm_rRNA_tRNA_snoRNA_miRNA_mRNA.fastq temp.txt


    total=`cat ../ready_files/"$A" | wc -l`
    rRNA=`cat "$A"_rRNA.fastq | wc -l`
    tRNA=`cat "$A"_tRNA.fastq | wc -l`
    snoRNA=`cat "$A"_snoRNA.fastq | wc -l`
    miRNA=`cat "$A"_miRNA.fastq | wc -l`
    mRNA=`cat "$A"_mRNA.fastq | wc -l`
    others=`cat "$A"_others.fastq | wc -l`
    unmappable=`cat "$A"_unmappable.fastq | wc -l`

    let total=total/4
    let rRNA=rRNA/4
    let tRNA=tRNA/4
    let snoRNA=snoRNA/4
    let miRNA=miRNA/4
    let mRNA=mRNA/4
    let others=others/4
    let unmappable=unmappable/4

    rm temp.txt
    rm "$A"_rm_rRNA.fastq
    rm "$A"_rm_rRNA_tRNA.fastq
    rm "$A"_rm_rRNA_tRNA_snoRNA.fastq
    rm "$A"_rm_rRNA_tRNA_snoRNA_miRNA.fastq
    rm "$A"_rm_rRNA_tRNA_snoRNA_miRNA_mRNA.fastq



printf "%s\t%d\t%d\t%d\t%d\t%d\t%d     \t%d\t%d\n" $A $total $rRNA $tRNA $snoRNA $miRNA $mRNA $others $unmappable >> ../analysis_results/small_RNA_profile.txt

done

cd ..




