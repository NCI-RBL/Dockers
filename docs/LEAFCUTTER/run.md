### Run STAR:

Align your samples in twopassMode for more accurate junction calling:

```bash
STAR \
  --runThreadN 64 \
  --runMode alignReads \
  --genomeDir 'STAR_hg38' \
  --sjdbGTFfile 'gencode.v44.annotation.gtf' \
  --outFileNamePrefix ${sample}_ \
  --outSAMunmapped Within \
  --outFilterType BySJout \
  --outFilterMultimapNmax 20 \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverLmax 0.04 \
  --alignIntronMin 20 \
  --alignIntronMax 50000 \ 
  --alignMatesGapMax 1000000 \
  --alignSJoverhangMin 5 \
  --alignSJDBoverhangMin 3 \
  --sjdbScore 1 \
  --readFilesIn ${sample}_R1.fq \
                ${sample}_R2.fq \
  --outFilterMatchNminOverLread 0.66 \
  --outSAMtype BAM Unsorted \
  --quantMode TranscriptomeSAM \
  --peOverlapNbasesMin 10 \
  --alignEndsProtrude 10 ConcordantPair \
  --twopassMode Basic \
  --outSAMstrandField intronMotif 

samtools sort -@ 64 -o ${sample}_Aligned.sorted.bam ${sample}_Aligned.out.bam
samtools index -@ 64 -M ${sample}_Aligned.sorted.bam 

```


### Run regtools


```bash
for bamfile in `ls *_Aligned.sorted.bam`; do
    echo Converting $bamfile to $bamfile.junc
    regtools junctions extract -s XS -a 8 -m 50 -M 500000 $bamfile -o $bamfile.junc
done
```

Repeat alignment for all samples. Leafcutter documentation recommends 4 samples per group mimimum.

### Download and Run leafcutter container

```bash
export SINGULARITY_CACHEDIR=$PWD/.singularity
singularity pull --dir ./ leafcutter.sif docker://wilfriedguiblet/leafcutter:v0.1
```

Create a juncfiles.txt. Example:
```
Sample1_Aligned.sorted.bam.junc
Sample2_Aligned.sorted.bam.junc
Sample3_Aligned.sorted.bam.junc
Sample4_Aligned.sorted.bam.junc
Sample5_Aligned.sorted.bam.junc
Sample6_Aligned.sorted.bam.junc
Sample7_Aligned.sorted.bam.junc
Sample8_Aligned.sorted.bam.junc
```

Download custom leafcutter:

git clone https://github.com/wilfriedguiblet/leafcutter.git


Enter the container
```bash
singularity shell -B $(pwd)/:/data2/,/data/RBL_NCI/:/data/RBL_NCI/ leafcutter.sif
```

Cluster juncfiles
```bash
python leafcutter/clustering/leafcutter_cluster_regtools.py -j juncfiles.txt -m 50 -o Experiment -l 500000
```

Prep annotation for leafcutter
```bash
leafcutter/leafviz/gtf2leafcutter.pl -o gencode.v44 gencode.v44.annotation.gtf
```


Create groups_file.txt. Example:
```
Sample1_Aligned.sorted.bam WT
Sample2_Aligned.sorted.bam WT
Sample3_Aligned.sorted.bam WT
Sample4_Aligned.sorted.bam WT
Sample5_Aligned.sorted.bam KO
Sample6_Aligned.sorted.bam KO
Sample7_Aligned.sorted.bam KO
Sample8_Aligned.sorted.bam KO
```

Run leafcutter:
```bash
Rscript leafcutter/scripts/leafcutter_ds.R --num_threads 64 Experiment_perind_numers.counts.gz groups_file.txt --min_samples_per_intron 4 --min_samples_per_group 4 -o Experiment
```

Create RData file for leafviz:
```bash
Rscript ../../leafcutter/leafviz/prepare_results.R --meta_data_file groups_file.txt \
  --code leafcutter Experiment_perind_numers.counts.gz \
  Experiment_cluster_significance.txt \
  Experiment_effect_sizes.txt \
  annotation_codes/gencode_v44/gencode.v44 \
  -o Experiment.RData
```

Exit container. Run leafviz locally (for instance with RStudio)
```r
options(shiny.host = "0.0.0.0")
options(shiny.port = 3838)
library(leafviz)
leafviz("~/Lab_Work/CCRRBL-5/leafviz/Experiment.RData")
```






