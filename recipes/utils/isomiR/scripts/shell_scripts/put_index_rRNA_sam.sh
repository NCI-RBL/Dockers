
ls -m1 *.sam | sed -e 's/\..*$//' > FF


cat FF | while read A
do
    echo "processing $A..."
    sh ~/Desktop/NGS/scripts/shell_scripts/put_index.sh ./$A

done

rm FF









