"""Complete supplementary comparison; original campaign stays immutable."""
import json
import statistics
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
MAIN = Path('/Users/theodosiou/Documents/Projects/test_barbac')
FINAL = MAIN / '.codex/publication-final-2026-09-12/source/benchmark/publication_final'
sys.path.insert(0, str(FINAL))
from common import load, save, sha
from analyze import paired_summary

METHODS = ['hamming','lv','shepherd_documented','starcode_sphere','starcode_mp','bartender']
METRICS = ['fn','fp','f1_percent','positive_truth_fn','positive_truth_f1_percent','incorrect_reads',
           'unassigned_reads','read_assignment_accuracy_percent','abundance_total_variation','zero_read_truth']


def main():
    receipt = load(HERE/'execution_validation.json')
    assert receipt['cells'] == 121
    original = load(FINAL/'timing_results.json')
    fresh = load(HERE/'results.json')
    retained = []
    for row in original:
        if row['method'] == 'shepherd':
            continue
        copy = dict(row)
        copy['reported_workflow_seconds'] = row['selected_workflow_seconds']
        copy['comparison_source'] = 'original campaign; accepted original or registered sleep-repair time'
        retained.append(copy)
    for row in fresh:
        copy = dict(row)
        copy['reported_workflow_seconds'] = row.get('workflow_seconds')
        copy['comparison_source'] = 'post hoc Shepherd configuration; later-session timing'
        retained.append(copy)
    data = pd.DataFrame(retained)
    assert len(data) == 726 and not data.duplicated(['condition','seed','method']).any()
    summaries = []
    for condition in ['random_mixed','anchored_mixed','milos']:
        for method in METHODS:
            rows = [r for r in retained if r['condition']==condition and r['method']==method]
            expected = 1 if condition=='milos' else 60
            assert len(rows) == expected
            ok = [r for r in rows if r['status']=='complete']
            complete = len(ok)==expected
            summary = dict(condition=condition,method=method,planned_n=expected,n_successful=len(ok),
                           complete_population=complete)
            for metric in METRICS:
                summary[metric] = statistics.mean(r[metric] for r in ok) if complete else None
                summary[metric+'_sd'] = statistics.stdev(r[metric] for r in ok) if complete and expected>1 else None
                if complete:
                    pandas_mean = data.loc[(data.condition==condition)&(data.method==method),metric].mean()
                    assert abs(summary[metric]-pandas_mean)<1e-10
            summary['workflow_seconds'] = statistics.median(r['reported_workflow_seconds'] for r in ok) if complete else None
            summaries.append(summary)
    contrasts = []
    for index, condition in enumerate(['random_mixed','anchored_mixed']):
        paired = data[(data.condition==condition)&data.method.isin(['lv','shepherd_documented'])]
        assert len(paired)==120
        if not (paired.status=='complete').all():
            contrasts.append(dict(condition=condition,status='incomplete',reason='Full 60-pair comparison unavailable'))
            continue
        values = paired.pivot(index='seed',columns='method',values='f1_percent')
        assert len(values)==60 and not values.isna().any().any()
        result = paired_summary(values.lv-values.shepherd_documented,2026091301+index)
        result['interpretation'] = 'Post hoc supplementary contrast; not retroactively confirmatory'
        contrasts.append(dict(condition=condition,status='complete',**result))
    save(HERE/'combined_results.json',retained)
    data.to_csv(HERE/'combined_results.csv',index=False)
    save(HERE/'summary.json',summaries)
    pd.DataFrame(summaries).to_csv(HERE/'benchmark_table.csv',index=False)
    save(HERE/'accuracy_contrasts.json',contrasts)
    pd.DataFrame(contrasts).to_csv(HERE/'accuracy_contrasts.csv',index=False)
    save(HERE/'analysis_validation.json',dict(rows=len(data),summaries=len(summaries),
        all_original_other_method_cells_reused=True,original_sources_sha256={str(FINAL/n):sha(FINAL/n)
          for n in ['results.json','timing_results.json','accuracy_contrasts.json','benchmark_table.pdf']},
        complete_shepherd_calls=len([r for r in fresh if r['status']=='complete']),
        summary_means_independently_reconciled=True,
        statistical_scope='Post hoc configuration sensitivity; no final-seed Barbac tuning'))
    print(json.dumps(contrasts,indent=2))


if __name__ == '__main__':
    main()
