#Aug 3rd, 2017
# No Parameters
# Convert all tsv files in tsv fold to a new fold named txt



echo "convert tsv to txt..."



mkdir txt

ls ./tsv/ -m1 *.* > samples

cat samples | while read one
do
	echo "processing $one..."
	~/Desktop/NGS/scripts/Tail_Trim_Analysis_v1.exe ./tsv/$one ./txt/"$one".txt
done


	
