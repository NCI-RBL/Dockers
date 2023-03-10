# sh adaptor_remove.sh min_length_passing_filter
# May 2017, changed the name from adaptor_remove.sh to the current one.
# deal with 8nt 3'adaptor sequences
# you need to input the minimal length for passing the filter.
# July 2021, updated the 3' adaptor sequence (8nt), the new exe file is now v7


echo "trimming adaptors..."



mkdir -p ready_files

ls -m1 *.fastq > samples

cat samples | while read one
do
	echo "processing $one..."
	/data/RBL_NCI/Gu/miRNA_NGS/scripts/adaptor_remove_read1_Chunmei $1 < $one > ./ready_files/"$one"_ready
done

echo `date` >> log.txt

