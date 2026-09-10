# Source records

The study is Jasinska et al. (2020), *Chromosomal barcoding of E. coli populations
reveals lineage diversity dynamics at high resolution*,
[doi:10.1038/s41559-020-1103-z](https://doi.org/10.1038/s41559-020-1103-z).

- `supplementary_tables.xlsx` is the publisher's
  [supplementary workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41559-020-1103-z/MediaObjects/41559_2020_1103_MOESM3_ESM.xlsx).
  The accompanying CSVs are extracted sheets. Table 1c and Table 4b are the
  selected constant-concentration experiment; Table 1b is retained for the
  initial experiment-selection audit and is not used in the final comparison.
- `PRJNA592529_runs.tsv` and `constant_aliases.tsv` preserve ENA run metadata
  for the selected constant-concentration experiment.
- `PRJNA592371_runs.tsv` and `increasing_aliases.tsv` preserve metadata for the
  increasing-concentration experiment considered during selection, not analyzed.
- The representative BioSample XML records explicitly link well identifiers to
  named treatments and replicates. `well_conditions.csv` is their parsed mapping;
  `prepare_sources.R` rebuilds the selected sample manifests from these records.
- The full paper was read from the publicly accessible
  [author-hosted copy](https://websites.umich.edu/~zhanglab/clubPaper/11_01_2022.pdf).
  The cached PDF and extracted text are local generated files, not repository
  artifacts. The cassette FASTA is transcribed from the Methods reference.

These records were retrieved on 10 September 2026. The sample manifests retain
ENA run accessions, URLs, byte counts and MD5 checksums. Processing receipts
retain the reference and source SHA-256 hashes and verified input MD5s.
`extraction_validation.json` records a synthetic indel/orientation check through
the actual mapping and barbac extraction path. `search_validation.json` records
the v13/v14 assignment equivalence check, not a controlled speed benchmark.
