# Tailer workflow

This tutorial describes how to use [Tailer](https://github.com/TimNicholsonShaw/tailer) on High-performance Cluster such as ***Biowulf*** or ***FRCE***, then [Tailer-analysis](https://github.com/TimNicholsonShaw/tailer-analysis) locally on RStudio.



## Install Tailer

You may want to run this into a dedicated environment.

```BASH
module load python
pip install jla-tailer
```


## Run Tailer

```BASH
Tailer -a annotation.gtf sample.bam
```


## Local Tailer-analysis

Download and extract [Tailer-analysis](https://hpc.nih.gov/~RBL_NCI/tailer-analysis-local.zip)

Open in RStudio. Click on Run App.



