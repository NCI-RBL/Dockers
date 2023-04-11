### Run pipeline:

Simply submit the pipeline as a slurm job:

```bash
sbatch isomiR.sh <ExperimentName> <Trimmer> <MinLength> <Consensus>
```

Trimmer options:

- Qiagen
- Illumina
- NEB
- GuLab


MinLength of trimming. Requires an integer.


Consensus: path to "motif-consensus.fa".
