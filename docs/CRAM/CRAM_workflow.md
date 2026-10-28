# cram.sh Tutorial

`cram.sh` converts alignment files in either direction based on their filename:

- `sample.bam` becomes `sample.cram`
- `sample.cram` becomes `sample.bam`

It checks files with `samtools quickcheck`, coordinate-sorts inputs whose header does not declare `SO:coordinate`, and writes an index by default.

## Cluster Setup

Download the script from [GitHub](https://github.com/NCI-RBL/Dockers/blob/main/workflows/CRAM/cram.sh).

The script requires Bash and `samtools` on `PATH`. On clusters that use environment modules, load samtools in the shell or job environment, for example:

```bash
module load samtools
```

Keep the script in a persistent location accessible from the login and compute nodes, then make it executable:

```bash
chmod +x /path/to/cram.sh
```

For convenient interactive use, add an alias to `~/.bashrc`, substituting the actual script path:

```bash
alias CRAM='/path/to/cram.sh'
```

Load the alias in the current shell with:

```bash
source ~/.bashrc
```

Aliases are generally not available in scheduler batch jobs. In a job script, call `/path/to/cram.sh` directly and load samtools there.

## References

Every conversion needs a reference FASTA and its samtools index (`.fai`). The script recognizes these built-in names:

| Name | FASTA path |
| --- | --- |
| `hg38` | `/fdb/genomebrowser/fasta/hg38/hg38.fa` |
| `mm10` | `/fdb/genomebrowser/fasta/mm10/mm10.fa` |

Use a built-in name with `-r`:

```bash
CRAM -r hg38 sample.bam
CRAM -r mm10 sample.cram
```

Or provide a FASTA path directly:

```bash
CRAM -r /path/to/reference.fa sample.bam
```

The FASTA must exist and be indexed. If needed, create the index with:

```bash
samtools faidx /path/to/reference.fa
```

A reference can also be set for repeated commands in the current shell:

```bash
export CRAM_REFERENCE=hg38
CRAM sample.bam
```

Use `CRAM_THREADS` and `CRAM_VERSION` to set the corresponding defaults:

```bash
export CRAM_THREADS=8
export CRAM_VERSION=3.1
```

`CRAM_VERSION` controls BAM-to-CRAM conversion. The default thread count uses `nproc` when available, falling back to 4.

## Common Workflows

Convert one or more BAMs to CRAM:

```bash
CRAM -r hg38 sample1.bam sample2.bam
```

Convert CRAM back to BAM:

```bash
CRAM -r hg38 sample1.cram
```

Convert all matching files in the current directory and put outputs in another directory:

```bash
CRAM -r hg38 -o converted *.bam
```

Use a specific thread count and verify input/output record counts:

```bash
CRAM -r mm10 -t 8 -c sample.bam
```

## Options

| Option | Meaning |
| --- | --- |
| `-r REF` | Built-in reference name (`hg38`, `mm10`) or indexed FASTA path. |
| `-t INT` | Threads per samtools call. |
| `-v STR` | CRAM version for BAM-to-CRAM conversion (`3.0` or `3.1`; default `3.1`). |
| `-o DIR` | Write outputs to this directory; it is created if needed. |
| `-f` | Replace existing output files. Without this option, existing outputs are skipped. |
| `-n` | Do not write an output index (`.crai` or `.bai`). |
| `-c` | Compare input and output record counts before completing. |
| `-d` | Delete the input and its recognized index files after successful conversion and count verification. Implies `-c`. Use with care. |
| `-h` | Show command-line help. |

Without `-o`, each output is placed next to its input. Inputs that do not end in `.bam` or `.cram` are skipped. An input whose header lacks `SO:coordinate` is sorted by coordinate before conversion; already coordinate-sorted inputs are not re-sorted. Temporary files are removed after successful conversion.

## Reference Caution

The script checks that the reference and its `.fai` exist, but that does not prove the reference is the correct assembly for an input. Samtools may detect a mismatch using reference checksums in CRAM metadata, but this is not guaranteed for every file. The `-c` option verifies record counts only; it does not verify that the reference sequence is correct.

## Troubleshooting

- **`samtools not found in PATH`**: load the samtools module or activate the environment that provides samtools.
- **`Reference not found`**: check the built-in mapping or pass the correct FASTA path with `-r`.
- **`Reference index missing`**: create the `.fai` using `samtools faidx`.
- **Output already exists**: choose another `-o` directory or use `-f` to replace it.
- **Sorting fails**: check available temporary disk space. `samtools sort` uses `$TMPDIR` when set, otherwise the output directory, for temporary sort data.

Run `CRAM -h` for the current built-in help.
