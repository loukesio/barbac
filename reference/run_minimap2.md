# Map merged reads to reference using minimap2 and process with samtools

Map merged reads to reference using minimap2 and process with samtools

## Usage

``` r
run_minimap2(
  merged_dir = "merged",
  reference,
  output_dir = file.path(merged_dir, "bam")
)
```

## Arguments

- merged_dir:

  Directory where merged fastq files are stored (default: "merged")

- reference:

  Path to reference FASTA file

- output_dir:

  Directory for BAM files (default: "merged/bam")

## Value

Vector of commands run
