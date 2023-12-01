### Run pipeline:

On biowulf,

Load the latest version of the pipeline:
```bash
export PATH=$PATH:/data/RBL_NCI/iCLIP/latest/
```

Run the pipeline with genomic coordinates:
```bash
sbatch iCLIP_latest.sh "your_directory"
```

After running with genomic coordinates, you can also run for transcriptomic coordinates:
```bash
sbatch iCLIP_Transcriptome_latest.sh "your_directory"
```

You can add the option `-resume` to resume a failed/interrupted run.
