"""Explicit post-start reporting addendum; retain failures, test only complete pairs."""
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.stats import t

from analyze import holm, paired_summary
from common import HERE, load, save, sha, verify_outputs
from run import check_freeze

COMPETITORS = ['shepherd', 'starcode_sphere', 'starcode_mp', 'bartender']
METRICS = ['fn', 'fp', 'f1_percent', 'positive_truth_f1_percent', 'incorrect_reads',
           'unassigned_reads', 'read_assignment_accuracy_percent',
           'abundance_total_variation', 'input_reads', 'zero_read_truth']


def summarize(cells, condition, method, expected):
    success = cells[cells.status == 'complete']
    row = dict(condition=condition, method=method, n=len(cells),
               n_successful=len(success), n_failed=len(cells)-len(success),
               planned_n=expected, complete_population=len(success)==expected)
    for metric in METRICS:
        row[metric] = float(success[metric].mean()) if len(success) else None
        row[metric+'_sd'] = float(success[metric].std(ddof=1)) if len(success)>1 else None
    row.update(workflow_seconds=float(success.workflow_seconds.median()) if len(success) else None,
        workflow_q25=float(success.workflow_seconds.quantile(.25)) if len(success) else None,
        workflow_q75=float(success.workflow_seconds.quantile(.75)) if len(success) else None)
    conditional = dict(row, scope='Conditional on successful execution; not a full-design comparison')
    if not row['complete_population']:
        for metric in METRICS:
            row[metric] = row[metric+'_sd'] = None
        for metric in ['workflow_seconds', 'workflow_q25', 'workflow_q75']:
            row[metric] = None
    return row, conditional


def compare(final, protocol, resamples=100000):
    n = protocol['n_independent_libraries_per_design']
    seeds = protocol['final_seeds']
    contrasts, speeds = [], []
    for condition in protocol['designs']:
        selected = final[(final.condition == condition) & (final.status == 'complete')]
        f1 = selected.pivot(index='seed', columns='method', values='f1_percent').reindex(index=seeds, columns=protocol['methods'])
        times = selected.pivot(index='seed', columns='method', values='workflow_seconds').reindex(index=seeds, columns=protocol['methods'])
        for competitor in COMPETITORS:
            complete = f1[['lv', competitor]].notna().all(axis=1)
            base = dict(condition=condition, competitor=competitor, planned_n=n, n=int(complete.sum()))
            if not complete.all():
                base.update(status='unavailable_missing_pairs',
                            missing_seeds=[int(s) for s in complete.index[~complete]])
                contrasts.append(dict(**base, mean_difference=None, sd_difference=None,
                    se_difference=None, mean_ci95_lower=None, mean_ci95_upper=None,
                    simultaneous_one_sided_lower=None, bootstrap_lower=None, one_sided_p=None,
                    wins=None, ties=None, losses=None, superiority_supported=False))
                speeds.append(dict(**base, geometric_time_ratio=None, ci95_lower=None, ci95_upper=None))
                continue
            assert np.isfinite(times[['lv', competitor]]).all().all()
            assert (times[['lv', competitor]] > 0).all().all()
            result = paired_summary(f1.lv-f1[competitor], 20260912401+len(contrasts),
                                    family_size=8, resamples=resamples)
            contrasts.append(dict(**{k:v for k,v in base.items() if k!='n'},
                                  status='complete', **result))
            logs = np.log(times.lv/times[competitor])
            half = t.ppf(.975, n-1)*logs.std(ddof=1)/np.sqrt(n)
            speeds.append(dict(**base, status='complete', geometric_time_ratio=float(np.exp(logs.mean())),
                ci95_lower=float(np.exp(logs.mean()-half)), ci95_upper=float(np.exp(logs.mean()+half)),
                inference_scope='Secondary descriptive paired log-time interval; no primary multiplicity guarantee'))
    assert len(contrasts) == 8
    adjusted = holm([r['one_sided_p'] if r['status']=='complete' else 1.0 for r in contrasts])
    for row, value in zip(contrasts, adjusted):
        row['holm_adjusted_p'] = float(value) if row['status']=='complete' else None
    return contrasts, speeds


def main():
    check_freeze()
    protocol = load(HERE/'final_protocol.json')
    validation = load(HERE/'execution_validation.json')
    rows = load(HERE/'results.json')
    data = pd.DataFrame(rows)
    n = protocol['n_independent_libraries_per_design']
    expected = {(s,c,m) for s in protocol['final_seeds'] for c in protocol['designs'] for m in protocol['methods']}
    final = data[data.scope == 'final']
    actual = {(int(r.seed),r.condition,r.method) for r in final.itertuples()}
    assert actual == expected and len(final)==len(expected)
    fixed = data[data.scope == 'published_reference']
    assert len(fixed)==6 and set(fixed.method)==set(protocol['methods']) and set(fixed.condition)=={'milos'}
    assert len(rows)==validation['cells']==n*12+6
    assert data.status.isin(['complete','failed']).all()
    assert int((data.status=='failed').sum())==validation['failed_cells']
    failures=[]
    for row in rows:
        key='milos' if row['scope']=='published_reference' else f"{int(row['seed'])}/{row['condition']}"
        dest=HERE/'generated/results'/key/row['method']
        assert load(dest/'result.json')==row
        verify_outputs(dest,row['output_sha256'])
        if row['status']=='failed':
            logs={name:(dest/name).read_text() for name in ['tool.log','worker.log'] if (dest/name).exists()}
            failures.append(dict(condition=row['condition'],seed=row['seed'],method=row['method'],
                error=row['error'],command=row['command'],logs=logs,
                preserved_sha256={name:sha(dest/name) for name in ['result.json','tool.log','worker.log','input.tsv'] if (dest/name).exists()}))
    contrasts,speeds=compare(final,protocol)
    summary,conditional=[],[]
    for condition in protocol['designs']+['milos']:
        for method in protocol['methods']:
            cells=data[(data.condition==condition)&(data.method==method)]
            expected_n=1 if condition=='milos' else n
            assert len(cells)==expected_n
            whole,successful=summarize(cells,condition,method,expected_n)
            summary.append(whole)
            conditional.append(successful)
    for name,value in [('accuracy_contrasts',contrasts),('speed_contrasts',speeds),
                       ('summary',summary),('successful_only_summary',conditional)]:
        save(HERE/(name+'.json'),value)
        pd.DataFrame(value).to_csv(HERE/(('benchmark_table' if name=='summary' else name)+'.csv'),index=False)
    data.to_csv(HERE/'all_results.csv',index=False)
    save(HERE/'failure_audit.json',failures)
    save(HERE/'analysis_validation.json',dict(primary_contrasts=8,
        available_contrasts=sum(r['status']=='complete' for r in contrasts),
        independent_libraries_per_design=n,rows=len(rows),failed_cells=len(failures),
        supported_contrasts=sum(r['superiority_supported'] for r in contrasts),
        all_final_cells_included=True,milo_inferentially_excluded=True,
        plan_unchanged=False,reporting_addendum='FAILURE_REPORTING.md',
        original_analysis_unchanged=True,complete_contrast_method_unchanged=True,
        family_size_retained=8,no_failed_scores_imputed=True,no_reduced_sample_inference=True,
        all_receipts_and_successful_output_hashes_verified=True,
        addendum_sha256=sha(HERE/'FAILURE_REPORTING.md'),driver_sha256=sha(HERE/'analyze_available.py'),
        completed_at=datetime.now(timezone.utc).isoformat()))
    print('RETAINED FAILURES',len(failures))
    print('ACCURACY CONTRASTS',[(r['condition'],r['competitor'],r['mean_difference'],r['status'],r['superiority_supported']) for r in contrasts])


if __name__=='__main__':
    main()
