# Changelog

## barbac 0.2.0

- Publish the validated v14 clustering engine on the main development
  line, including exact LV candidate partitions, abundance-bound
  pruning, support-based tie ordering and the optional Poisson indel
  model.
- Preserve observed barcode lengths with flank-based BAM extraction,
  with optional read identifiers and extraction diagnostics.
- Support R1-only and mixed single/paired sample tables in
  [`run_cli_pipeline()`](https://loukesio.github.io/barbac/reference/run_cli_pipeline.md).
  PEAR is required only for overlapping paired reads. Required command
  failures stop the run, existing output directories must be empty, and
  the returned sample table links original labels to their indexed BAMs.
- Use the minimap2 short-read preset and retain primary alignments in
  the pipeline; mapping statistics count primary reads or merged
  molecules.
- Bundle all 32 LTC palettes for direct use in barbac plots.
- Add barbac Studio: a local Shiny app for count tables or FASTQ input,
  native clustering, lineage plots, statistics, downloads and Quarto
  report exports.
- Include synthetic video walkthroughs, a concise front-page guide, and
  linked scientific workflows and validation evidence. Historical
  benchmark receipts retain their original engine versions and
  measurements.

## barbac 0.1.0

- Initial DNA barcode clustering, extraction, plotting and CLI workflow.
