# Package index

## Core clustering

Abundance-ranked greedy Levenshtein clustering with a bit-parallel C++
kernel.

- [`super_cluster2()`](https://loukesio.github.io/barbac/reference/super_cluster2.md)
  : Fast Centroid-Based Sequence Clustering

## Time-series visualisation

Stacked-area plots of lineage frequencies over time.

- [`barbac_ts_area()`](https://loukesio.github.io/barbac/reference/barbac_ts_area.md)
  : Stacked-area plot of barcode lineage frequencies over time
- [`theme_barbac()`](https://loukesio.github.io/barbac/reference/theme_barbac.md)
  : A publication-ready ggplot2 theme for barbac plots

## Barcode extraction

- [`barbac_xtr()`](https://loukesio.github.io/barbac/reference/barbac_xtr.md)
  : Extract Barcodes from BAM File
- [`barbac_xtr.stats()`](https://loukesio.github.io/barbac/reference/barbac_xtr.stats.md)
  : Generate Summary Statistics and Plots for Extracted Barcodes

## End-to-end pipeline orchestrator

- [`run_cli_pipeline()`](https://loukesio.github.io/barbac/reference/run_cli_pipeline.md)
  : Run Full Barbac CLI Pipeline

## Individual pipeline steps

Wrappers around FastQC, PEAR, minimap2 and samtools.

- [`run_fastqc()`](https://loukesio.github.io/barbac/reference/run_fastqc.md)
  : Run FastQC on all FASTQ files listed in a samples.csv file
- [`run_multiqc()`](https://loukesio.github.io/barbac/reference/run_multiqc.md)
  : Run MultiQC to summarize FastQC reports
- [`run_pear_merge()`](https://loukesio.github.io/barbac/reference/run_pear_merge.md)
  : Merge paired-end reads using PEAR
- [`run_minimap2()`](https://loukesio.github.io/barbac/reference/run_minimap2.md)
  : Map merged reads to reference using minimap2 and process with
  samtools
- [`summarise_bam_stats()`](https://loukesio.github.io/barbac/reference/summarise_bam_stats.md)
  : Summarise mapped and unmapped read counts from sorted BAM files

## Environment management

Provision and activate a conda environment with the CLI tools.

- [`configure_environment()`](https://loukesio.github.io/barbac/reference/configure_environment.md)
  : Configure Conda Environment for barbac
- [`use_barbac_env()`](https://loukesio.github.io/barbac/reference/use_barbac_env.md)
  : Use barbac conda environment in the current R session
- [`check_barbac_tools()`](https://loukesio.github.io/barbac/reference/check_barbac_tools.md)
  : Check availability of tools in barbac conda environment

## Package overview

- [`barbac-package`](https://loukesio.github.io/barbac/reference/barbac-package.md)
  [`barbac`](https://loukesio.github.io/barbac/reference/barbac-package.md)
  : barbac: End-to-end DNA barcode lineage analysis
