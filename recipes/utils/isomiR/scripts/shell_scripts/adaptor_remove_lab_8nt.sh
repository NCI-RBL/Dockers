# sh adaptor_remove.sh min_length_passing_filter
# May 2017, changed the name from adaptor_remove.sh to the current one.
# deal with 8nt 3'adaptor sequences
# you need to input the minimal length for passing the filter.


echo "trimming adaptors..."



mkdir ready_files

ls -m1 *.* > samples

cat samples | while read one
do
	echo "processing $one..."
	~/Desktop/NGS/scripts/adaptor_remove_read1_v5.exe $1 < $one > ./ready_files/"$one"_ready
done

echo `date` >> log.txt
	
