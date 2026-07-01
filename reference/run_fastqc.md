# Run FastQC on all FASTQ files listed in a samples.csv file

This function reads a \`samples.csv\` file from a given path, extracts
\`R1\` and optionally \`R2\` FASTQ paths, and runs FastQC on all unique
files. The output is saved under \`results/fastQC/\` in the same
directory as the CSV.

## Usage

``` r
run_fastqc(sample_csv)
```

## Arguments

- sample_csv:

  Path to the \`samples.csv\` file.

## Value

A character vector of FastQC system commands executed.
