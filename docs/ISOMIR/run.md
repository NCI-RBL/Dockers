### Run pipeline:

Simply submit the pipeline as a slurm job:

```bash
sbatch isomiR.sh <ExperimentName> <Trimmer> <MinLength> <Consensus> <Index>
```

Trimmer options:

- Qiagen
- Illumina
- NEB
- GuLab


MinLength of trimming. Requires an integer.


Consensus: path to "motif-consensus.fa".


Currently available indexes:

- hg38
- mm39
