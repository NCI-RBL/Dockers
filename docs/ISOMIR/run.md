### Run pipeline:

Simply submit the pipeline as a slurm job:

```bash
sbatch isomiR.sh <ExperimentName> "<Options>"
```

Example run:

```bash
sbatch isomiR.sh Name "--trimmer Qiagen --MinLen 18 --consensus motif-consensus.fa --index hg38"
```

List of options.

* --trimmer : Trimmer used. Default is Qiagen.
    * Qiagen
    * Illumina
    * NEB
    * GuLab
    * None

* --MinLen: Minimum length for trimming. Requires an integer. Default is 18.

* --consensus : alternative consensus file.

* --index : Index for profiling. Currently available:
    * hg38
    * mm39

* --sRNAprofiling : Yes or No. No skips this step.

* --QuagmiR : Yes or No. No skips this step.
