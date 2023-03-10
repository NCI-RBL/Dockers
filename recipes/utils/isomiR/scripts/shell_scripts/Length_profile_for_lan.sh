# sh Length_profile_for_lan.sh

ls -m1 *.fastq > samples

cat samples | while read one
do
	~/Desktop/NGS/scripts/Len_profile_fastq.exe < $one > "$one".txt
done

	
