# Extract Barcodes from BAM File

This function extracts barcode sequences from a BAM file based on the
specified reference name and position range, and saves the result as a
clean CSV file.

## Usage

``` r
barbac_xtr(
  bam_file,
  ref_name = "Reference_barcodes",
  start_pos = 54,
  end_pos = 78,
  output_file = NULL,
  min_count = 1,
  verbose = TRUE
)
```

## Arguments

- bam_file:

  Character string. Path to the BAM file from which to extract barcodes.

- ref_name:

  Character string. Name of the reference to be used. Default is
  "Reference_barcodes".

- start_pos:

  Numeric. Start position of the barcode in the reference. Default is
  54.

- end_pos:

  Numeric. End position of the barcode in the reference. Default is 78.

- output_file:

  Character string or NULL. Path for output CSV. If NULL, creates
  filename based on input BAM file. Default is NULL.

- min_count:

  Numeric. Minimum count threshold to include a barcode. Default is 1.

- verbose:

  Logical. Print progress messages. Default is TRUE.

## Value

A character string with the path to the generated CSV file.

## Details

This function:

- Reads aligned sequences from a BAM file at specified genomic
  coordinates

- Extracts barcode sequences from the specified region

- Counts occurrence of each unique barcode

- Filters by minimum count threshold

- Saves results to CSV with columns: barcode, counts, barcode_length

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract barcodes from default region
barbac_xtr("sample1_sorted.bam")

# Extract from custom region
barbac_xtr(
  bam_file = "sample1_sorted.bam",
  ref_name = "chr1",
  start_pos = 100,
  end_pos = 125
)

# Specify output file and filter low counts
barbac_xtr(
  bam_file = "sample1_sorted.bam",
  output_file = "sample1_filtered_barcodes.csv",
  min_count = 5
)
} # }
```
