# Complete-sample mapping pilot, 9 September 2026

**Historical pilot:** the complete eight-sample BAM extraction and comparison
have since finished. See [current results](results/README.md). The measurements
and planning estimates below describe the earlier mapping-only check.

The selected eight samples are practical to start locally. The first complete
sample, **SRR22757105 / hBFA1-YPD-R1-Time8**, ran through the existing
`run_cli_pipeline()` in **248.38 seconds (4.14 minutes)**, including R startup,
FastQC, PEAR, minimap2, BAM generation/indexing/statistics and MultiQC.
Download, barcode extraction and clustering are excluded. This is one run,
not a repeated runtime benchmark.

The machine has 16 GiB RAM and eight physical/logical CPU cores. The run used
the existing `barbac_env` and default pipeline thread settings; PEAR used one
thread. No new environment was installed. Full-process peak memory was not
recorded: macOS denied the timing wrapper's resource query after printing its
elapsed/CPU times. The R pipeline completed, and its BAM and QC outputs were
verified. The wrapper's nonzero status must not be presented as a measured
memory result. [The receipt](mapping_pilot.json) retains this distinction.

| Outcome | Count |
|---|---:|
| Input read pairs, complete files verified against ENA sizes/MD5s | 1,367,175 |
| Merged reads produced by PEAR | 1,365,286 |
| Pairs not merged | 1,889 |
| Primary mapped merged reads | 1,364,530 |
| Unmapped merged reads | 756 |
| Secondary or supplementary alignments | 0 |

**99.945% of merged reads mapped**, equivalent to **99.807% of input pairs**
producing a mapped merged read. These are mapping retention percentages, not
barcode accuracy, recovery of true lineages, or F1. The BAM passes
`samtools quickcheck`; FastQC reports for both mates and the MultiQC HTML exist.

The second-resolution pipeline log attributes about 13 seconds to FastQC,
172 to PEAR, 27 to mapping/BAM processing, one to BAM statistics and 30 to
MultiQC. These sum to 243 seconds; the outer elapsed measurement also includes
R startup. The selected eight runs contain 16,511,755 read pairs. Simple scaling
by read count suggests approximately **50 minutes for sequential FASTQ-to-BAM
processing**, assuming similar cost per pair. That is a planning estimate;
the complete analysis, extraction and method comparison have not been timed.

## Reference and extraction status

[The candidate FASTA](reference/chen2023_masked_amplicon.fasta) is a 167-base
engineered amplicon with two 26-base N blocks. It was reconstructed from modal
fixed segments in the first 50,000 read pairs and verified against the reverse
complement of the shared middle segment in the other mate. It uses no published
barcode identities. [Provenance](reference/provenance.json) records the evidence.
It is **not an author-deposited FASTA**, and the full sample includes the prefix
used for reconstruction; this is a mapping feasibility check, not independent
biological validation.

Nominal coordinates are BC2 50–75 and BC1 110–135 (one-based inclusive). BC1
is reverse-complemented relative to the author's read-oriented barcode key.
The older project reference describes a different barcode construct.

Before comparing methods, the BAM extraction needs to preserve observed
insertions/deletions, both component identities, original read-pair quality
filtering and the first-passing-observation UMI rule. The existing fixed-width
`barbac_xtr()` call alone does not implement that protocol. This pilot stops at
BAM generation. The earlier direct extractor remains a publication-parser
validation baseline; its SLURM stage is not yet the requested final BAM workflow.

## Reproduce with the existing barbac environment

Use the normal barbac installation and `configure_environment()` setup from
the main README. Then download the checked first sample into a fresh directory:

```bash
python3 benchmark/time_series_chen2023/run_sample.py download \
  --index 0 --work-dir /absolute/path/to/chen_work
Rscript benchmark/time_series_chen2023/run_mapping_pilot.R \
  /absolute/path/to/chen_work
```

The R helper calls `use_barbac_env()` and `run_cli_pipeline()`, verifies required
outputs and BAM integrity, and saves session information and an elapsed-time
receipt. It rejects an existing output directory. It does not run barcode
extraction or clustering. Reconstruct the reference from the saved 50,000-pair
pilot with `python3 benchmark/time_series_chen2023/build_mapping_reference.py`.
