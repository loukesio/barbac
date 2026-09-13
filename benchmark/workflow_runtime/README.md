# Real-data FASTQ-to-plot runtime

**33 minutes 18 seconds** from local FASTQ files through saved PDF and PNG plots
on an Apple M1 with 16 GiB RAM. The run used **10,754,210 reads**, all 12 technical
runs from the first three A3 timepoints of the Jasińska 2020 series. Native LV
clustering took **343.514 seconds (5 minutes 44 seconds)**. Tables were available
in approximately ten minutes; drawing every lineage in both formats took a
further 23 minutes 18 seconds.

| Stage | Elapsed seconds |
|---|---:|
| Package loading | 0.725 |
| Combine technical runs | 0.447 |
| FASTQ → BAM and QC | 203.783 |
| Barcode extraction | 45.638 |
| Pool counts | 1.566 |
| Native LV clustering | 343.514 |
| Cluster statistics | 0.008 |
| Sample assignment and table exports | 3.725 |
| All-lineage PDF + PNG exports | 1,398.070 |
| **Outer R process, total** | **1,997.977** |

The outer time includes startup and finishing overhead, so it is slightly larger
than the sum of named stages. Command-level FastQC, mapping and MultiQC timings
are nested within the FASTQ-to-BAM stage; do not add them again.

The workflow mapped 8,591,566 reads and extracted **8,578,419 barcode reads**,
representing **437,009 distinct sequences** and **144,537 inferred lineages**.
Every extracted count was preserved through clustering and sample assignment.
All lineages are plotted; no abundance cutoff or remainder category is used.
The 15-base nominal barcodes and dense set of nearby sequences make this a
different search workload from Milo. Input size alone does not predict runtime.

This is one descriptive application run, with ordinary desktop/browser activity,
not a replicated speed comparison. No other agent-launched analysis jobs ran
concurrently. No system sleep/wake event overlapped the run. Minimap2 retained
its three-thread default; clustering and BLAS used one thread. Downloads,
installation and checksum verification are excluded; input concatenation is
included. The raw-read example does not apply the study-specific Q10 filter or
UMI deduplication. It is not a reproduction of the separate filtered full-study
analysis.

## Reproduce

Obtain the public FASTQs with the Jasińska application download script and verify
the [study manifest](../time_series_jasinska2020/samples.tsv). Prepare the current
isolated release, then choose a new output directory:

```sh
Rscript app/run.R --prepare-only
python3 benchmark/workflow_runtime/run.py \
  --input-root benchmark/time_series_jasinska2020 \
  --output benchmark/workflow_runtime/generated/my-run
```

The launcher uses Studio's active release library; `--library` can select an
explicit installation. The exact measured code revision is `ecc26af3e1a4c983e0cc3110fd4e455247634565`,
recorded before the run. Later launcher changes protect installation folders;
the timed R pipeline, extraction, clustering and plotting code are unchanged.

[Recorded selection and timing scope](protocol.md) ·
[Stage times](results/stages.csv) · [Commands](results/command_timings.csv) ·
[Full receipt](results/receipt.json) · [Count summary](results/result.json) ·
[Timing verification](results/timing_validation.json)

The raw measured plots are retained locally in `generated/run-01`. A separate
presentation correction shortens a clipped vertical axis label for the manuscript;
its formatting time is outside this recorded workflow. Memberships, frequencies
and plotted data are unchanged. Large generated inputs and result tables remain
local; the scripts and compact receipts are versioned.
