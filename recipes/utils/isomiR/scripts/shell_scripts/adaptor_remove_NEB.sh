# sh adaptor_remove_illumina.sh
# new scripts to handle results of libarary prepared by illumina kit
# used an old program - adaptor_remove_illumina_small_rna_kit.exe
# 2016-May-10th changed to trim NEB adaptors

echo "trimming adaptors... the minimal length of reads passing filter is 10nt"



mkdir ready_files

ls -m1 *.* > samples

cat samples | while read one
do
	echo "processing $one..."
	~/Desktop/NGS/scripts/adaptor_remove_NEB_v1.exe < $one > ./ready_files/"$one"_ready
done

echo `date` >> log.txt
	
