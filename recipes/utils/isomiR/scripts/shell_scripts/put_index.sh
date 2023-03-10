~/Desktop/NGS/tools/samtools-1.8/samtools view -b -S $1.sam > $1.bam
~/Desktop/NGS/tools/samtools-1.8/samtools sort $1.bam -o $1.sorted.bam
~/Desktop/NGS/tools/samtools-1.8/samtools index $1.sorted.bam
