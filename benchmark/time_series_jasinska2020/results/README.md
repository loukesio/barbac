# Reconstructed E. coli time series

The [explained Quarto report](../r_report/report.html) presents the results.
The [analysis instructions](../README.md) describe inputs, processing and SLURM.

The 75 longitudinal samples contain 350,934,958 reads, of which 305,326,997
(87.00%) yield extracted barcodes. The three shared initial samples contain
13,332,953 reads and yield 10,533,102 extracted barcode reads. Count the shared
baseline once when reporting unique sequencing input, even though it is reused
within each population's retrospective clustering.

Across the 120 published dominant-barcode entries, the median population-level
Spearman correlation is 0.9985. Mean absolute frequency differences are 0.548
percentage points using all input reads and 0.066 points using extracted barcode
reads. These are two explicitly different comparisons, not accuracy estimates;
the publication's denominator discrepancy is explained below.

At passage 30, median Shannon effective diversity is 19.5 lineages under
chloramphenicol and 97.5 without antibiotic. Individual treated replicates range
from 11.2 to 130.2, versus 46.3 to 113.8 in controls. The direction of the median
contrast is consistent with the study's low-chloramphenicol result, but the
replicates overlap and this is not a formal estimate of the diversity-loss rate.

| File | Contents |
|---|---|
| `A3`, `B3`, `C3` `_barcode_counts.csv.gz` | Complete chloramphenicol replicate 1, 2 and 3 count matrices |
| `A1`, `B1`, `C1` `_barcode_counts.csv.gz` | Complete untreated replicate 1, 2 and 3 count matrices |
| `extraction_summary.csv` | Read accounting and processing times for 78 samples, including the three initial samples |
| `diversity.csv` | Full-count diversity and native `cluster_stats()` summaries at every population/timepoint |
| `publication_sample_comparison.csv` | Matched sample summaries from barbac and Supplementary Table 1c |
| `publication_top_barcodes.csv` | Original all-input-normalized exports for the 120 published population/barcode entries |
| `publication_frequency_comparison.csv` | The same entries with both all-input and extracted-read normalizations |
| `publication_frequency_summary.csv` | Per-population rank correlation and mean/maximum absolute differences under both normalizations |
| `publication_denominator_audit.csv` | Consistency check between published extraction fractions and summed final frequencies |
| `clustering_times.csv` | Original development-build clustering timings (`-O0` in this run); archived, not used for reported speed |
| `release_clustering_times.csv` | Optimized release clustering elapsed and CPU seconds, including support ordering |
| `release_validation.json` | Compiler command, binary/source hashes and exact full-population equivalence of release/development results |
| `metadata_schema_migration.json` | Audit of the checksum-map repair; original cache signatures verified and every measured timing preserved |
| `processing_receipts.json` | Input checksums, source/reference hashes, counts and timings for every processed sample |
| `processing_validation.json` | Independent reconciliation of all 312 FastQC read totals and 78 receipts |
| `provenance.json` | Recorded settings, build marker and analysis-source hashes |
| `validation.json` | Independent count-matrix reconciliation and hashes of the result files |

The checked-in release timings were measured on a 16 GB Apple silicon Mac with
eight physical CPUs, using three concurrent release-clustering processes. Other
local analysis tasks were also active. They measure the `super_cluster2()` call,
including support ordering, and exclude input preparation, extraction, count
matrix export and report generation. They are application measurements, not a
controlled comparison with another method. The validation receipt records the
compiler command, R platform, binary hash and exact input/output checks.
Elapsed clustering times range from 1,238.237 to 1,517.868 seconds
(20.6–25.3 minutes) per population; active CPU times range from 1,220.242 to
1,497.021 seconds. All six release reruns reproduce every original membership
and centroid count exactly.

Count-matrix rows are observed centroid barcode identifiers; columns are sampled
passages, including the shared `passage_0`. Counts are sequencing reads, not cells
or UMI-deduplicated molecules. Zero is an observed zero count. Separate clustering
within each pooled population can assign shared baseline reads differently, so
baseline diversity can differ between reconstructions without representing
different starting biological populations.

For positive lineage counts `c`, the main diversity panels use `q = c / sum(c)`:
richness is the number of positive entries, Shannon effective diversity is
`exp(-sum(q * log(q)))`, and inverse dominant-lineage share is `1 / max(q)`.
These are the `assigned_*` fields. The unprefixed `shannon_effective` and
`inverse_dominance` exports instead substitute `p = c / input_reads`; when
extraction is incomplete, `p` does not sum to one. They record the alternative
calculation and must not be interpreted as standard effective-number indices.

The author frequency denominator remains unresolved. Low CMP r1's final top-20
frequencies sum to 86.94%, exceeding Table 1c's reported 80.24% extraction yield.
The two tables therefore cannot both use all input reads for that sample.
Both normalizations remain explicit in the report. Spearman correlation is
invariant to positive within-population rescaling; absolute differences are not.
Neither statistic is an accuracy percentage or a substitute for known truth.

Reference comparisons use exact centroid identifiers and retain zero calls.
An absent exact identifier can also result from assigning its reads to another
centroid; it is not automatically evidence that a biological lineage is absent.
The published list is a fixed set of abundant barcodes, not a truth set for all
rare lineages. Complete author count matrices and processing code would be
needed to reconcile every trajectory and the normalization discrepancy.
