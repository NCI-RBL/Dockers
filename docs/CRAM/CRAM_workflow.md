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

BAM-to-CRAM conversion needs a reference FASTA and its samtools index (`.fai`). CRAM-to-BAM conversion can use embedded reference slices or samtools' configured reference lookup without an explicit FASTA. The script recognizes these built-in names:

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

When a FASTA is supplied, it must exist and be indexed. If needed, create the index with:

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

Embed reference slices in CRAM output:

```bash
CRAM -r hg38 -e sample.bam
```

Convert CRAM back to BAM:

```bash
CRAM -r hg38 sample1.cram
```

For a CRAM with embedded reference slices, or a reference available through samtools' lookup, an explicit FASTA is optional. If `CRAM_REFERENCE` was set earlier, unset it to use automatic lookup:

```bash
unset CRAM_REFERENCE
CRAM sample1.cram
```

If `CRAM_REFERENCE` is set, the script uses it as the reference unless `-r` overrides it.

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
| `-r REF` | Built-in reference name (`hg38`, `mm10`) or indexed FASTA path. Required for BAM input unless `CRAM_REFERENCE` is set; optional for CRAM input when samtools can resolve the reference. |
| `-t INT` | Threads per samtools call. |
| `-v STR` | CRAM version for BAM-to-CRAM conversion (`3.0` or `3.1`; default `3.1`). |
| `-o DIR` | Write outputs to this directory; it is created if needed. |
| `-f` | Replace existing output files. Without this option, existing outputs are skipped. |
| `-n` | Do not write an output index (`.crai` or `.bai`). |
| `-e` | Embed reference slices in CRAM output; applies to BAM inputs only. |
| `-c` | Compare input and output record counts before completing. |
| `-d` | Delete the input and its recognized index files after successful conversion and count verification. Implies `-c`. Use with care. |
| `-h` | Show command-line help. |

Without `-o`, each output is placed next to its input. Inputs that do not end in `.bam` or `.cram` are skipped. An input whose header lacks `SO:coordinate` is sorted by coordinate before conversion; already coordinate-sorted inputs are not re-sorted. Temporary files are removed after successful conversion.

## Reference Validation

When a reference is supplied with `-r` or `CRAM_REFERENCE`, the script compares each input's `@SQ` sequence dictionary with it before conversion. It fails if a contig is missing, its length differs, or its `M5` checksum differs. `M5` is an MD5 checksum of a contig's sequence; it is checked only when present in the input header.

If an input lacks `M5` values, the script warns that it checked contig names and lengths only. Assemblies with identical names and lengths but different bases cannot be distinguished from that header. The `-c` option compares record counts only and does not establish reference identity.

## Embedded References

Without `-e`, BAM-to-CRAM conversion does not embed reference slices by default. With `-e`, samtools stores the reference portions needed by the CRAM slices, not the entire reference genome. These CRAMs can be decoded without an external FASTA.

When decoding CRAM without `-r` or `CRAM_REFERENCE`, samtools resolves references as it decodes each slice, using embedded slices or its configured reference lookup. The script does not perform a separate whole-file check for an embedded reference before decoding. If any required reference cannot be found, conversion fails and the partial output is removed. Supply the correct indexed FASTA with `-r` if needed.

Embedding normally increases CRAM size, but the overhead depends on the coverage, alignment layout, reference and compression settings. Measure it on a representative BAM using separate output directories:

```bash
CRAM -r hg38 -o external sample.bam
CRAM -r hg38 -e -o embedded sample.bam
wc -c external/sample.cram embedded/sample.cram
```

The percentage size increase is `(embedded bytes - external bytes) / external bytes * 100`. Compare the CRAM files themselves, not disk usage rounded to filesystem blocks.

## Troubleshooting

- **`samtools not found in PATH`**: load the samtools module or activate the environment that provides samtools.
- **`Reference not found`**: check the built-in mapping or pass the correct FASTA path with `-r`.
- **`Reference index missing`**: create the `.fai` using `samtools faidx`.
- **`reference check failed`**: use the same reference assembly used for the input; check contig names, lengths and any reported `M5` mismatch.
- **`Unable to fetch reference` or CRAM decoding fails**: provide the matching indexed FASTA with `-r`, or configure samtools reference lookup. A CRAM may need external reference data even when some slices have embedded references.
- **Output already exists**: choose another `-o` directory or use `-f` to replace it.
- **Sorting fails**: check available temporary disk space. `samtools sort` uses `$TMPDIR` when set, otherwise the output directory, for temporary sort data.

Run `CRAM -h` for the current built-in help.
