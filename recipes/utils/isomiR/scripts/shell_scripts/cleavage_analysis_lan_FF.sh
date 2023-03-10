
echo $1 >> $1.txt


echo	>> $1.txt
echo "All reads, threshold 1..." >> $1.txt
echo	>> $1.txt

~/Desktop/NGS/scripts/cleavage_site_lan_analysis.exe ~/Desktop/NGS/index/Lan/FF_luc_target.fa $1 all no 1 off >> $1.txt
echo   >> $1.txt

for i in {6..25}
do
	~/Desktop/NGS/scripts/cleavage_site_lan_analysis.exe ~/Desktop/NGS/index/Lan/FF_luc_target.fa $1 $i no 1 off >> $1.txt
done


echo	>> $1.txt
echo "All reads, threshold 200..." >> $1.txt
echo	>> $1.txt

~/Desktop/NGS/scripts/cleavage_site_lan_analysis.exe ~/Desktop/NGS/index/Lan/FF_luc_target.fa $1 all no 200 off >> $1.txt
echo   >> $1.txt

for i in {6..25}
do
	~/Desktop/NGS/scripts/cleavage_site_lan_analysis.exe ~/Desktop/NGS/index/Lan/FF_luc_target.fa $1 $i no 200 off >> $1.txt
done

echo	>> $1.txt
echo	>> $1.txt
echo "Unique reads, threshold 1..."	>> $1.txt
echo	>> $1.txt

~/Desktop/NGS/scripts/cleavage_site_lan_analysis.exe ~/Desktop/NGS/index/Lan/FF_luc_target.fa $1 all yes 1 off >> $1.txt
echo   >> $1.txt

for i in {6..25}
do
	~/Desktop/NGS/scripts/cleavage_site_lan_analysis.exe ~/Desktop/NGS/index/Lan/FF_luc_target.fa $1 $i yes 1 off >> $1.txt
done

echo	>> $1.txt
echo	>> $1.txt
echo "Unique reads, threshold 200..."	>> $1.txt
echo	>> $1.txt

~/Desktop/NGS/scripts/cleavage_site_lan_analysis.exe ~/Desktop/NGS/index/Lan/FF_luc_target.fa $1 all yes 200 off >> $1.txt
echo   >> $1.txt

for i in {6..25}
do
	~/Desktop/NGS/scripts/cleavage_site_lan_analysis.exe ~/Desktop/NGS/index/Lan/FF_luc_target.fa $1 $i yes 200 off >> $1.txt
done