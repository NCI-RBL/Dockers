
	adaptor_remove_read1_v4.exe $1 < $one > ./ready_files/"$one"_ready

mkdir ready_files

echo "trimming adaptors..."


	echo "processing $one..."
cat samples | while read one
