
# you need to input the minimal length for passing the filter.


echo "trimming adaptors..."
echo "you need to input the minimal length for passing the filter!"


mkdir ready_files

ls -m1 *.* > samples

cat samples | while read one
do
	echo "processing $one..."
	~/Desktop/NGS/scripts/adaptor_remove_Qiagen_v3 $1 < $one > ./ready_files/"$one"_ready
done

echo `date` >> log.txt
	
