"""Predeclared paired-library accuracy contrasts; Milo is descriptive only."""
import json

import numpy as np
import pandas as pd
from scipy.stats import t
from common import HERE, load, save


def holm(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    order = np.argsort(pvalues, kind='stable')
    adjusted = np.empty(len(pvalues))
    adjusted[order] = np.minimum(1, np.maximum.accumulate(pvalues[order] * (len(pvalues) - np.arange(len(pvalues)))))
    return adjusted


def paired_summary(difference, seed, family_size=8, resamples=100000):
    difference = np.asarray(difference, dtype=float)
    n = len(difference)
    if n < 2 or not np.isfinite(difference).all():
        raise ValueError('Need finite differences from independent paired libraries')
    mean, sd = float(difference.mean()), float(difference.std(ddof=1))
    se = sd / np.sqrt(n)
    if se == 0:
        p = 0.0 if mean > 0 else 1.0
    else:
        p = float(t.sf(mean / se, n - 1))
    lower = float(mean - t.ppf(1 - .05 / family_size, n - 1) * se)
    rng = np.random.default_rng(seed)
    # Resample whole libraries, preserving each within-library method pair.
    means = np.empty(resamples)
    for start in range(0, resamples, 5000):
        size = min(5000, resamples-start)
        means[start:start+size] = difference[rng.integers(0, n, (size, n))].mean(axis=1)
    bootstrap_lower = float(np.quantile(means, .05 / family_size))
    return dict(n=n, mean_difference=mean, sd_difference=sd, se_difference=float(se),
        mean_ci95_lower=float(mean-t.ppf(.975, n-1)*se), mean_ci95_upper=float(mean+t.ppf(.975, n-1)*se),
        simultaneous_one_sided_lower=lower, bootstrap_lower=bootstrap_lower,
        one_sided_p=p, wins=int((difference > 1e-10).sum()), ties=int((np.abs(difference) <= 1e-10).sum()),
        losses=int((difference < -1e-10).sum()), superiority_supported=bool(lower > 0 and bootstrap_lower > 0))


def main():
    protocol = load(HERE / 'final_protocol.json')
    validation = load(HERE / 'execution_validation.json')
    if validation['failed_cells']:
        raise RuntimeError('Report retained tool failures before attempting complete paired inference')
    rows = load(HERE / 'results.json')
    data = pd.DataFrame(rows)
    n = protocol['n_independent_libraries_per_design']
    final = data[data['scope'] == 'final']
    assert len(final) == n * 2 * 6
    assert not final.duplicated(['seed', 'condition', 'method']).any()
    contrasts, speeds = [], []
    for condition in protocol['designs']:
        selected = final[final.condition == condition]
        f1 = selected.pivot(index='seed', columns='method', values='f1_percent')
        times = selected.pivot(index='seed', columns='method', values='workflow_seconds')
        assert len(f1) == n and not f1.isna().any().any()
        for competitor in ['shepherd', 'starcode_sphere', 'starcode_mp', 'bartender']:
            delta = f1.lv - f1[competitor]
            result = paired_summary(delta, 20260912401 + len(contrasts))
            contrasts.append(dict(condition=condition, competitor=competitor, **result))
            logs = np.log(times.lv / times[competitor])
            half = t.ppf(.975, n-1) * logs.std(ddof=1) / np.sqrt(n)
            speeds.append(dict(condition=condition, competitor=competitor, n=n,
                geometric_time_ratio=float(np.exp(logs.mean())), ci95_lower=float(np.exp(logs.mean()-half)),
                ci95_upper=float(np.exp(logs.mean()+half)), inference_scope='Secondary descriptive paired log-time interval; no primary multiplicity guarantee'))
    for row, adjusted in zip(contrasts, holm([r['one_sided_p'] for r in contrasts])):
        row['holm_adjusted_p'] = float(adjusted)
    save(HERE / 'accuracy_contrasts.json', contrasts)
    pd.DataFrame(contrasts).to_csv(HERE / 'accuracy_contrasts.csv', index=False)
    save(HERE / 'speed_contrasts.json', speeds)
    pd.DataFrame(speeds).to_csv(HERE / 'speed_contrasts.csv', index=False)
    summary = []
    for condition in protocol['designs'] + ['milos']:
        for method in protocol['methods']:
            x = data[(data.condition == condition) & (data.method == method)]
            assert len(x) == (1 if condition == 'milos' else n)
            row = dict(condition=condition, method=method, n=len(x))
            for metric in ['fn', 'fp', 'f1_percent', 'positive_truth_f1_percent', 'incorrect_reads', 'unassigned_reads',
                           'read_assignment_accuracy_percent', 'abundance_total_variation', 'input_reads', 'zero_read_truth']:
                row[metric] = float(x[metric].mean())
                row[metric+'_sd'] = float(x[metric].std(ddof=1)) if len(x)>1 else None
            row.update(workflow_seconds=float(x.workflow_seconds.median()),
                workflow_q25=float(x.workflow_seconds.quantile(.25)), workflow_q75=float(x.workflow_seconds.quantile(.75)))
            summary.append(row)
    save(HERE / 'summary.json', summary)
    pd.DataFrame(summary).to_csv(HERE / 'benchmark_table.csv', index=False)
    metrics = [k for k in data.columns if k not in ['command', 'output_sha256', 'signature']]
    data[metrics].to_csv(HERE / 'all_results.csv', index=False)
    save(HERE / 'analysis_validation.json', dict(primary_contrasts=len(contrasts), independent_libraries_per_design=n,
        rows=len(rows), supported_contrasts=sum(r['superiority_supported'] for r in contrasts),
        all_final_cells_included=True, milo_inferentially_excluded=True, plan_unchanged=True))
    print('ACCURACY CONTRASTS', [(r['condition'], r['competitor'], round(r['mean_difference'],6), r['superiority_supported']) for r in contrasts])


if __name__ == '__main__':
    main()
