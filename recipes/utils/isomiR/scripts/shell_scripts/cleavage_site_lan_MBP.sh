mkdir cleavage_site_5P
mkdir cleavage_site_3P

ls -m1 *.sam | sed -e 's/\..*$//' > MBP
cat MBP | while read A
do
    echo "processing $A..."

    ~/Desktop/NGS/scripts/cleavage_site_lan.exe 10 ~/Desktop/NGS/index/Lan/MBP_target.fa "$A".sam ./cleavage_site_5P/"$A"_cleavage_site_5P.txt
    ~/Desktop/NGS/scripts/cleavage_site_3P_lan.exe 10 ~/Desktop/NGS/index/Lan/MBP_target.fa "$A".sam ./cleavage_site_3P/"$A"_cleavage_site_3P.txt
done

cd ..







