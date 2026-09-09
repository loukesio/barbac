### Barcode time-series application

We analyzed eight barcode-sequencing samples from the hBFA1 YPD fitness assay
of Chen, Johnson, Hérissant et al. (2023), covering generations 8, 16, 24 and 40
in two biological assay replicates. Of 16,511,755 paired reads, the common
mapping, variable-length barcode extraction and UMI-deduplication workflow
retained 16,051,344 molecules. We pooled each barcode component across the
eight samples for retrospective clustering at distance three and reconstructed
sample-specific paired-barcode counts using a consistent membership map.
Barbac LV with support ordering and the optional Poisson indel model achieved
median Spearman correlation 0.999968 with published counts over
2,314 published pair identifiers; the corresponding Hamming value was
0.999929. Median combined-component clustering workflow times
were 14.65 s and 11.97 s,
respectively, over three serial repeats. These are agreement measurements
against a differently processed experimental reference, not ground-truth
accuracy estimates. Starcode showed slightly higher reference agreement,
while Hamming barbac had the lowest median clustering time (Table 3).
The median proportion of LV molecules assigned to exact published pair IDs
was 83.09%; two abundant absent pairs account for
69.47% of the pooled unmatched mass,
and both have BC2 at least seven edits from every published BC2. Their absence
cannot be attributed to clustering error from this comparison alone.

**Table 3. Publication agreement and clustering workflow time for the Chen 2023 time series.**
Spearman correlations use all 2,314 published pair IDs, filling absent calls with
zero. Molecule coverage is the percentage of method-input molecules assigned to
those exact IDs; both agreement columns are medians across eight samples.
Times are medians of three serial repeats, summing the BC2 and BC1 workflows
including startup and required exports, excluding shared FASTQ preprocessing.
Shepherd uses the supplied error rate 0.005 because automatic estimation failed
on BC1; the value matches the pre-existing barbac configuration and was not
selected using publication agreement. Unassigned totals span all eight samples.

| Method | Median Spearman ρ | Median molecule mass on published pair IDs | Clustering, median seconds | Unassigned molecules, all samples |
|---|---:|---:|---:|---:|
| barbac LV + Poisson | 0.999968 | 83.09% | 14.65 | 0 |
| barbac Hamming | 0.999929 | 82.90% | 11.97 | 0 |
| Shepherd (e=0.005) | 0.930542 | 80.30% | 27.02 | 189,126 |
| Starcode sphere | 0.999991 | 83.13% | 15.74 | 0 |
| Starcode message passing | 0.999990 | 83.14% | 14.84 | 0 |
| Bartender | 0.999931 | 82.93% | 58.48 | 0 |
