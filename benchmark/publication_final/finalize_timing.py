"""Select registered repair observations explicitly; preserve original timing exports."""
import math
import statistics

import pandas as pd
from scipy.stats import t

from common import HERE, load, save, sha, verify_outputs
from timing_repair import key


def main():
    protocol=load(HERE/'final_protocol.json')
    repair_protocol=load(HERE/'timing_repair_protocol.json')
    repair_validation=load(HERE/'timing_repair_validation.json')
    repairs={r['key']:r for r in load(HERE/'timing_repair_results.json')}
    registered={r['key'] for r in repair_protocol['cells']}
    assert set(repairs)==registered and repair_validation['cells']==len(registered)==28
    originals=load(HERE/'results.json')
    assert len(originals)==726
    rows=[]
    for original in originals:
        row=dict(original)
        row['original_workflow_seconds']=original.get('workflow_seconds')
        row['selected_workflow_seconds']=original.get('workflow_seconds')
        row['timing_status']='original_valid' if original['status']=='complete' else 'unavailable_tool_failure'
        if key(original) in repairs:
            repair=repairs[key(original)]
            dest=HERE/'generated/timing_repair'/key(original)
            assert load(dest/'result.json')==repair
            verify_outputs(dest,repair['output_sha256'])
            assert sha(HERE/'generated/results'/key(original)/'result.json')==repair['original_receipt_sha256']
            if repair['status']=='complete':
                assert repair['accuracy_mapping_byte_equal']
                for name in ['centroids.csv','members.csv']:
                    assert repair['output_sha256'][name]==original['output_sha256'][name]
                row['timing_status']='repair_valid'
                row['selected_workflow_seconds']=repair['workflow_seconds']
            else:
                row['timing_status']='unavailable_repair_failure'
                row['selected_workflow_seconds']=None
            row['timing_repair_receipt']=str(dest/'result.json')
        rows.append(row)
    summary=load(HERE/'summary.json')
    conditional=load(HERE/'successful_only_summary.json')
    time_summary=[]
    for row,conditional_row in zip(summary,conditional):
        cells=[r for r in rows if r['condition']==row['condition'] and r['method']==row['method']]
        times=[r['selected_workflow_seconds'] for r in cells if r['selected_workflow_seconds'] is not None]
        full=len(times)==row['planned_n']
        out=dict(condition=row['condition'],method=row['method'],planned_n=row['planned_n'],
            valid_time_n=len(times),repaired_time_n=sum(r['timing_status']=='repair_valid' for r in cells),
            timing_complete_population=full,workflow_seconds=statistics.median(times) if full else None,
            workflow_q25=float(pd.Series(times).quantile(.25)) if full else None,
            workflow_q75=float(pd.Series(times).quantile(.75)) if full else None)
        row.update({k:v for k,v in out.items() if k not in ['condition','method','planned_n']})
        conditional_row['workflow_seconds']=statistics.median(times) if times else None
        conditional_row['valid_time_n']=len(times)
        time_summary.append(out)
    contrasts=[]
    for condition in protocol['designs']:
        for competitor in ['shepherd','starcode_sphere','starcode_mp','bartender']:
            left={r['seed']:r['selected_workflow_seconds'] for r in rows if r['condition']==condition and r['method']=='lv'}
            right={r['seed']:r['selected_workflow_seconds'] for r in rows if r['condition']==condition and r['method']==competitor}
            assert left.keys()==right.keys() and len(left)==60
            available=[s for s in protocol['final_seeds'] if left[s] is not None and right[s] is not None]
            base=dict(condition=condition,competitor=competitor,n=len(available),planned_n=60,
                inference_scope='Secondary descriptive paired log-time interval; registered sleep-affected observations repaired')
            if len(available)!=60:
                contrasts.append(dict(**base,status='unavailable_missing_pairs',geometric_time_ratio=None,ci95_lower=None,ci95_upper=None))
                continue
            logs=[math.log(left[s]/right[s]) for s in available]
            mean=statistics.mean(logs);half=t.ppf(.975,59)*statistics.stdev(logs)/math.sqrt(60)
            contrasts.append(dict(**base,status='complete',geometric_time_ratio=math.exp(mean),
                ci95_lower=math.exp(mean-half),ci95_upper=math.exp(mean+half)))
    save(HERE/'timing_results.json',rows)
    save(HERE/'timing_summary.json',time_summary)
    save(HERE/'publication_summary.json',summary)
    save(HERE/'publication_successful_only_summary.json',conditional)
    save(HERE/'timing_contrasts.json',contrasts)
    pd.DataFrame(summary).to_csv(HERE/'benchmark_table.csv',index=False)
    pd.DataFrame(rows).to_csv(HERE/'publication_results.csv',index=False)
    pd.DataFrame(contrasts).to_csv(HERE/'timing_contrasts.csv',index=False)
    save(HERE/'timing_selection_validation.json',dict(original_cells=len(originals),registered_repairs=28,
        accepted_repairs=sum(r['timing_status']=='repair_valid' for r in rows),
        unavailable_repairs=sum(r['timing_status']=='unavailable_repair_failure' for r in rows),
        original_scores_and_measurements_retained=True,selection_uses_registered_repair_not_minimum_time=True))
    print('TIMING SELECTION VERIFIED',len(rows),'cells')


if __name__=='__main__':main()
