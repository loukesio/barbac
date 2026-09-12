"""Truth-based scoring, independent of clustering programs and simulation RNG."""
import numpy as np
import pandas as pd


def validate_data(inputs, truth, labels):
    assert list(inputs) == ['barcode', 'counts']
    assert list(truth) == ['barcode', 'true_count']
    assert list(labels) == ['member', 'true_barcode', 'read_count']
    assert inputs.barcode.is_unique and truth.barcode.is_unique
    assert not labels.duplicated(['member', 'true_barcode']).any()
    for data, column, minimum in [(inputs, 'counts', 1), (truth, 'true_count', 0), (labels, 'read_count', 1)]:
        assert not data.isna().any().any()
        assert ((data[column] >= minimum) & (data[column] == np.floor(data[column]))).all()
    assert inputs.barcode.str.fullmatch('[ACGT]+').all()
    assert truth.barcode.str.fullmatch('[ACGT]+').all()
    assert labels.true_barcode.isin(truth.barcode).all()
    pd.testing.assert_series_equal(labels.groupby('member').read_count.sum().sort_index(),
                                   inputs.set_index('barcode').counts.sort_index(), check_names=False)
    realized = labels.groupby('true_barcode').read_count.sum().reindex(truth.barcode, fill_value=0)
    difference = realized.to_numpy() - truth.true_count.to_numpy()
    assert int(inputs.counts.sum()) == int(truth.true_count.sum())
    return dict(parent_count_mismatches=int(np.count_nonzero(difference)),
                parent_count_absolute_difference=int(np.abs(difference).sum()),
                multi_parent_observed_sequences=int((labels.groupby('member').true_barcode.nunique() > 1).sum()))


def score(centroids, members, inputs, truth, labels):
    assert centroids.central_barcode.is_unique and members.member.is_unique
    assert not centroids.isna().any().any() and not members.isna().any().any()
    assert (centroids.sum_counts > 0).all()
    assert members.member.isin(inputs.barcode).all()
    assert members.central_barcode.isin(centroids.central_barcode).all()
    counts = inputs.set_index('barcode').counts
    assert np.array_equal(members.member_count.to_numpy(), members.member.map(counts).to_numpy())
    mapped_counts = members.groupby('central_barcode').member_count.sum().sort_index()
    pd.testing.assert_series_equal(mapped_counts, centroids.set_index('central_barcode').sum_counts.sort_index(),
                                   check_dtype=False, check_names=False)
    actual = set(truth.barcode)
    positive = set(truth.loc[truth.true_count > 0, 'barcode'])
    found = set(centroids.central_barcode)
    tp, fn, fp = len(found & actual), len(actual - found), len(found - actual)
    positive_tp, positive_fn, positive_fp = len(found & positive), len(positive - found), len(found - positive)
    origins = labels.merge(members[['member', 'central_barcode']], on='member', how='left', validate='many_to_one')
    correct = int(origins.loc[origins.true_barcode == origins.central_barcode, 'read_count'].sum())
    unassigned = int(origins.loc[origins.central_barcode.isna(), 'read_count'].sum())
    total = int(inputs.counts.sum())
    assert int(centroids.sum_counts.sum()) + unassigned == total
    difference = centroids.set_index('central_barcode').sum_counts.subtract(truth.set_index('barcode').true_count, fill_value=0)
    tv = (float(difference.abs().sum()) + unassigned) / (2 * total)
    assert 0 <= tv <= 1 + 1e-12
    return dict(tp=tp, fn=fn, fp=fp, f1_percent=200 * tp / (2 * tp + fn + fp),
                positive_truth_fn=positive_fn, positive_truth_fp=positive_fp,
                positive_truth_f1_percent=200 * positive_tp / (2 * positive_tp + positive_fn + positive_fp),
                correct_reads=correct, incorrect_reads=total - correct - unassigned,
                unassigned_reads=unassigned, read_assignment_accuracy_percent=100 * correct / total,
                abundance_total_variation=tv, true_barcodes=len(actual), zero_read_truth=len(actual - positive),
                inferred_barcodes=len(found), input_reads=total, input_sequences=len(inputs))
