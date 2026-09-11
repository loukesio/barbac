-- Independent report-side metric calculation from observed benchmark counts.
-- reference_runs is loaded from results.csv by build_report.py.
SELECT method, process_seconds, core_seconds, input_preparation_seconds,
       pipeline_seconds, build_id, centroids, unique_centroids,
       tp, fn, fp, observed_fn, absent_truth_recovered, ws,
       1.0 * tp / (tp + fp) AS precision,
       1.0 * tp / (tp + fn) AS recall,
       2.0 * tp / (2 * tp + fn + fp) AS f1,
       positive_truth_fn, positive_truth_fp,
       2.0 * (99591 - positive_truth_fn) /
           (2 * (99591 - positive_truth_fn) + positive_truth_fn + positive_truth_fp)
           AS positive_truth_f1,
       correct_reads, misassigned_reads, unassigned_reads,
       1.0 * correct_reads / input_reads AS read_assignment_accuracy,
       output_reads, consistent_parent_reads, consistent_parent_accuracy,
       input_reads, mapping_counts_reconcile
FROM reference_runs
ORDER BY pipeline_seconds ASC;
