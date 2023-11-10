### Run pipeline:

Simply submit the pipeline as a slurm job:

```bash
sbatch RiboFootPrint.sh
```

If your previous run was not completed and you wish to resume, resubmit the slurm job as:

```bash
sbatch RiboFootPrint.sh -resume
```


Running this pipeline will download a Docker container in a folder named "work". To remove the container, delete the folder.

```bash
rm -r work
```
