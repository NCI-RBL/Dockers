# Nov 18th 2016

# step 1. remove adaptor
# step 2. map to reference genome
# step 3. generate hit counts for each sgRNA

# input is raw data from miseq in the format of fastq.
# output is 1. ready file 2. map file 3. hit counts.


echo "removing adaptors"

mkdir ready_files
mkdir no_adaptor

ls -m1 *.* > samples

cat samples | while read one
do
    echo ""
    echo ""
    echo "processing $one..."
	~/Desktop/NGS/scripts/adaptor_remove_CRISPR_library_v3_activation.exe ./ready_files/"$one"_ready ./no_adaptor/"$one"_noadaptor < $one
done

echo `date` >> log.txt
rm samples

echo ""
echo ""
echo "Start processing all _ready files in ./read_files/ folder..."

mkdir analysis_files
mkdir analysis_results

cd ready_files

ls -m1 *_ready > samples

cd ../analysis_files

cat ../ready_files/samples | while read A
do
echo ""
echo ""
echo "processing $A..."
echo "mapping to sgRNA library allowing 1 mismatch (v 1 mode)..."

~/Desktop/NGS/tools/bowtie-1.1.1/bowtie -v 1 --best --norc --sam ~/Desktop/NGS/index/CRISPRa_SAM_V2-sgRNA  --un "$A"_unmapped.fastq  --al  "$A"_mapped.fastq  ../ready_files/"$A" "$A".sam

echo ""
echo "generating counts by sgRNA..."

~/Desktop/NGS/scripts/analyze_CRISPR_count.exe ~/Desktop/NGS/index/CRISPRa_SAM_V2-sgRNA_list.txt "$A".sam ../analysis_results/"$A"_count_by_sgRNA.txt

done

cd ..
