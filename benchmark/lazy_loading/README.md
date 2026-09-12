# Faster startup with unchanged v14 clustering

The candidate defers loading BAM-processing namespaces until extraction is used. The extraction calls already use explicit package-qualified functions, so only two redundant namespace imports and their roxygen declarations were removed. DESCRIPTION dependencies, public function signatures, clustering rules and native C++ source are unchanged.

## Paired full-workflow LV measurements

Three alternating baseline/candidate pairs per input, using the exact frozen Python worker and identical timing boundaries. Values below are medians; ranges and all paired observations are saved separately.

| Input | v14 seconds | Candidate seconds | Time reduction |
|---|---:|---:|---:|
| random_mixed | 4.923 | 1.914 | 61.1% |
| anchored_mixed | 5.146 | 2.041 | 60.3% |
| milos | 38.310 | 35.032 | 8.6% |

These timings establish the local baseline/candidate improvement. They do not replace frozen competitor measurements or establish a new fastest-tool ranking. The complete workflow includes worker/package startup, input handling, clustering and exports. The native clustering implementation is unchanged.

## Correctness and startup

- All 18 paired workflow outputs reproduce frozen v14 centroids, memberships and counts byte-for-byte.
- Ten additional compatibility runs cover both Hamming and LV across all five frozen inputs, including the substitution-only controls.
- All 526 assertions across 44 package tests passed, with zero failures, errors, skips or warnings.
- In five alternating fresh-startup pairs, median process time fell from 3.729 to 0.614 seconds. GenomicAlignments stayed unloaded until needed by extraction.
- Main, its released library, all competitor outputs and the reserved final seeds remain untouched.

One initial baseline pilot completed clustering but the timing logger failed while building its result record because a command field was duplicated. Its protocol and outputs are retained as workflow_protocol_attempt1.json and generated/paired. The reporting bug was fixed before the complete 18-run campaign; no pilot timing was selected into the reported results.

## Reproduction

Build this checkout into a separate R library, then run validate.R, measure.py and paired_workflows.py from this directory tree. The benchmark scripts record their exact local source/library paths. They expect empty output directories for fresh timings and preserve existing measurements. See protocol.json, startup_results.json, compatibility.json, workflow_protocol.json, workflow_observations.json, workflow_summary.json and both validation receipts.

The speed candidate is ready for integration review independently of the count-learning research candidate. It does not add a FASTQ, base-quality or FastQC requirement to super_cluster2().
