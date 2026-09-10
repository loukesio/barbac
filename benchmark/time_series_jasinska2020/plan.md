# E. coli application and explanatory Quarto report

## Scope agreed with user

Use R, existing barbac extraction/statistics/plotting functions, Quarto HTML,
gt summary tables and A–D diagnostic captions. This explicit surface overrides
the reporting skills' generic JavaScript artifact packaging. Audience: technical
methods, with sufficient experimental introduction for a reader new to the paper.
No competitor runs in this application; existing benchmark evidence stays separate.

Selected before processing: constant drug experiment PRJNA592529; chloramphenicol
1 microgram/mL and untreated control; three deposited biological replicates per
condition; all available timepoints and four sequencing subsamples per timepoint.
Shared initial population: three samples, four subsamples each, pooled for baseline.
Selection is based on a simple fixed treatment and complete deposited series, not
on achieved agreement. Each passage is approximately six generations in the source.

## Evidence and comparison

- Author paper: DOI 10.1038/s41559-020-1103-z; 288-base masked cassette, barcode
  reference positions 11–25 (15 bases), single-end sequencing.
- Sample annotations explicitly map wells to treatment/replicate. Subsamples
  are partitions of a sequencing sample, not extra biological replicates.
- Supplementary Table 1c: input reads, quality/extraction fractions and barcode
  richness by sample. Compare only exact matching treatment/replicate/passage.
- Supplementary Table 4b: published top barcode sequences, average and final
  frequencies. Full per-barcode per-timepoint author counts are not provided in
  these tables. Do not claim a complete trajectory-by-trajectory replication.
- Match the paper's minimum Phred 10 across a read. No UMI rule is described;
  counts are reads, not deduplicated molecules.
- Reads begin at/near the variable region. Minimap2 maps the downstream constant
  sequence; barbac_xtr selects that mapped anchor (32–40) and captures the upstream
  query barcode with ^([ACGT]{10,20})TATCTCGGTAG. Ns in the published reference
  are placeholders; reference alignment must not replace observed barcode bases.
  This explicit extraction rule differs from the author's custom alignment rule.
- Primary barbac: LV, distance 3, support order, ratio 20, assumed error rate .005,
  Poisson indel option, design scoring off. Pool within each population over time
  with shared initial data, following the retrospective trajectory design.
- Initial interpretation of the Methods: assigned lineage reads / all sample
  reads. Validation found that Low CMP r1's published final top-20 frequencies
  sum to 86.94%, exceeding its reported 80.24% extraction fraction. The author
  denominator therefore remains unresolved; both all-input and conditional
  extracted-read comparisons are retained, with explicitly labelled denominators.
  Calculate all diversity metrics from complete counts, never from a top-N plot.

## Report structure and chart contracts

Title and summary; study explanation and experimental timeline before detailed
findings; selected samples/definitions; extraction diagnostics; trajectories and
diversity; publication comparison; methods, limitations and next questions. This
reorders the technical-report skill's definitions before findings to address the
user's explicitly stated confusion. Audit files and commands stay in README.

- Experiment schematic: two conditions, three replicates, actual sampling times.
- Extraction figure: A length, B abundance (explain log10), C base-composition
  entropy in bits, D native gt length-summary table; distinct sequences count once.
- Sequence entropy is within one barcode. Author Extended Data 1 uses entropy
  across barcodes at each position (natural logarithm). Population diversity uses
  lineage abundances. Explain all three, do not present them as interchangeable.
- Lineage trajectories: actual observed points, common barcode identity/color,
  explicit zeros/missingness and denominator; selected before evaluating agreement.
- Population diversity: richness, exponential Shannon, inverse maximum frequency;
  actual timepoints, all three replicates visible, sampling-depth caveat.
- Publication comparison: per-sample richness and final top-barcode frequency
  comparisons with identity line and all underlying numeric values available.
- Tables: native gt for compact summary; searchable DT for detailed count lookup.

## Validation

ENA checksums/sizes; read counts vs manifest and author tables; quality/mapping/
extraction accounting; synthetic boundary extraction; complete member assignments
and count conservation; known diversity examples; zero handling; source hashes;
package tests for optional stats details/default compatibility; rendered offline
HTML inspection with interactive controls. No claim of known-truth accuracy or
antibiotic-resistance phenotype from barcode frequencies alone.
