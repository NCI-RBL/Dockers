/Users/gus6/Desktop/NGS/tools/samtools-1.2/samtools index $1.sorted.bam
/Users/gus6/Desktop/NGS/tools/samtools-1.2/samtools view -b -S $1.sam > $1.bam
/Users/gus6/Desktop/NGS/tools/samtools-1.2/samtools sort $1.bam $1.sorted
