# sh adaptor_remove_illumina.sh
# new scripts to handle results of libarary prepared by illumina kit
# used an old program - adaptor_remove_illumina_small_rna_kit.exe
# 2016-May-10th changed to trim NEB adaptors

#2016-Sep-6th: replace the exe file with trim_read1_mis_v4_5nt to handle library prepared with our own protocol.

echo "trimming adaptors... the minimal length of reads passing filter is 7nt"



mkdir ready_files

ls -m1 *.* > samples

cat samples | while read one
do
	echo "processing $one..."
	~/Desktop/NGS/scripts/trim_read1_mis_v4_5nt.exe 7 < $one > ./ready_files/"$one"_ready
done

echo `date` >> log.txt
	