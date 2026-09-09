# Comparison after complete extraction

The real-data question is whether each method reconstructs reproducible barcode
pair counts and trajectories, and how much compute it takes. The published
counts are a reference analysis, **not a known true set**. Differences from that
analysis must not be called false positives or false negatives.

1. Confirm all eight extraction jobs finish, input hashes match the manifest,
   molecule accounting balances, and inline indices/lengths match the study.
   Check UMI conflicts and quality losses before interpreting rare lineages.
2. Freeze these shared extracted inputs and parameters before comparing methods.
   Use distance three as in the existing benchmarks. Compare current barbac
   Hamming, LV with the explicitly labeled Poisson option, Shepherd, both
   Starcode modes, and Bartender. Record exact versions, full workflow elapsed
   time and peak resident memory on the same allocated hardware. Repeat timing
   measurements serially; separate extraction costs from clustering costs.
3. For a retrospective trajectory comparison, pool component counts across the
   eight selected samples separately for BC2 and BC1. For each method, use its
   resulting member-to-centroid maps to reconstruct pairs and assign each
   sample's original molecule counts. Retain paired identities throughout.
   Report unresolved molecules. This pooled analysis uses later time points;
   it does not establish prospective performance on unseen future samples.
4. Align resulting `BC2_BC1` identifiers with the author's
   `hBFA1_all_freqs_tidy.csv`, restricted to YPD. Validate orientation first.
   Compare counts and abundance-stratified agreement; count unmatched identifiers
   on both sides and the molecule mass they represent. These are agreement
   measures, not TP/FN/FP labels. Inspect high-count disagreements individually.
5. Report barcode-count Spearman correlations per sample and trajectory plots
   for selected abundant and rare shared pairs. For normalized comparisons,
   declare each denominator. Give shared-set conditional frequencies alongside
   the fraction of each sample's retained molecules represented by that set.
   Never silently normalize one method over all molecules and another over
   only matched centroids. The author's `Clipped_Log10(Freq)` contains clipping;
   it cannot be treated as an exact invertible frequency measurement.
6. Document the publication's pooled correction, >10-count filters, lane
   intersection, chimera removal and GC/timepoint filters. Differences in
   processing scope can explain differences in final counts. Do not optimize
   barbac parameters against the reference and then present that agreement as
   an independent accuracy evaluation.
7. Only after count-level agreement is understood, compare growth trends or
   fitness estimates using the same eligible time points and neutral-lineage
   normalization as the publication. Preserve genuine late/rare lineages when
   joining samples; do not restrict everything to the earliest observed set.

Deliverables: one sample × method accuracy-agreement/runtime table, matched and
unmatched molecule summaries, a trajectory figure, resource-use measurements,
and manuscript text distinguishing simulation truth from experimental agreement.
No claim that barbac is universally best follows from high extraction retention
or agreement with one published processing pipeline.

Execution adjustment: Shepherd's automatic error-rate estimator failed for the
low-diversity BC1 input. Use its documented `-e 0.005` option consistently on
both components, matching the pre-existing barbac configured rate. This is a
supplied assumption, not an estimated experimental rate or tuning against the
publication's results. Retain automatic-attempt logs separately; only successful
fixed-rate runs enter the three-repeat timing comparison.
