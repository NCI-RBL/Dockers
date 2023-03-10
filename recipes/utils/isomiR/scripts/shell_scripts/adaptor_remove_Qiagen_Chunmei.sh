
# you need to input the minimal length for passing the filter.


echo "trimming adaptors..."
echo "you need to input the minimal length for passing the filter!"

mkdir -p ready_files
/mnt/rnabl-work/Guiblet/CCRRBL10/NGS/scripts/adaptor_remove_Qiagen_v3 $1 < $2 > ./ready_files/"$2"_ready
	
