# Complete application timing, registered before execution

This single application measurement answers how long the released Barbac workflow
takes from local FASTQ files to saved lineage plots on this machine. It is separate
from the frozen comparison of clustering methods; competitors are not rerun.

Use the first three chronological samples (passages 2, 4 and 6; generations 12,
24 and 36) of well A3, low chloramphenicol replicate 1, in the archived Jasińska
2020 manifest. Include all four technical sequencing runs per sample without
subsampling. Verify ENA MD5 hashes before timing. Selection is chronological,
independent of measured runtime or clustering results.

Time one fresh R process using the isolated, optimized 0.2.1 release library.
Include package loading, concatenation of technical-run gzip files, FastQC,
the public `run_cli_pipeline()` with its unchanged default mapping settings,
MultiQC, `barbac_xtr()` with the archived anchored 10–20 nt extraction pattern,
count pooling, `super_cluster2()` with the publication LV settings, cluster
statistics, sample assignment, data exports, and `barbac_ts_area()` PDF and PNG
exports. Plot every inferred lineage with zero-filled missing counts, with no
abundance cutoff. Preserve raw counts throughout. Record per-stage and outer
process elapsed time, inputs, source revision, hashes, hardware and tool versions.

The generic raw-read example does not apply the publication-specific Q10 filter
used by the separate full Jasińska reanalysis, or deduplicate UMIs. It must not be
described as reproducing that filtered study analysis. Download, input verification
and package/environment installation are outside the measurement. Input staging
is inside it. CLI command timings are nested within the FASTQ-to-BAM stage and
must not be added to that stage a second time.

Keep this machine awake, run samples sequentially, limit BLAS/native clustering to
one thread, and launch no competing computational jobs during measurement.
The CLI retains minimap2's three-thread default and the other tools' defaults.
Record sleep assertions and sleep/wake logs. Retain every attempt: a failed or
sleep-interrupted attempt is labelled as such; any repair is explained and never
selected because it was faster. One run describes this workload and environment,
not a replicated speed comparison or a universal performance guarantee.
