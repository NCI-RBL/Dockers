# sh adaptor_remove_illumina.sh
# new scripts to handle results of libarary prepared by illumina kit
# used an old program - adaptor_remove_illumina_small_rna_kit.exe


echo "trimming adaptors... the minimal length of reads passing filter is 10nt"



mkdir -p ready_files

ls -m1 *.* > samples

cat samples | while read one
do
	echo "processing $one..."
	/home/malagobadans2/NGS/scripts/adaptor_remove_Ryan < $one > ./ready_files/"$one"_ready
done

echo `date` >> log.txt
	
