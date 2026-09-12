"""Independent arithmetic and base-R checks after all timed work has finished."""
import csv
from datetime import datetime
import math
from pathlib import Path
import statistics
import subprocess

from common import HERE, REFERENCE, load, save, sha
from run import check_freeze
from timing_repair import key, START, END


def close(a,b,absolute=1e-10):
    assert math.isclose(a,b,rel_tol=1e-10,abs_tol=absolute),(a,b)


def unique_csv(path,key,value,cast=str):
    data={}
    with path.open() as stream:
        for row in csv.DictReader(stream):
            assert row[key] not in data
            data[row[key]]=cast(row[value])
    return data


def check_mapping(source,dest,row):
    inputs=unique_csv(source/'input.csv','barcode','counts',int)
    truth=unique_csv(source/'truth.csv','barcode','true_count',int)
    centroids=unique_csv(dest/'centroids.csv','central_barcode','sum_counts',int)
    members=unique_csv(dest/'members.csv','member','central_barcode')
    total=sum(inputs.values())
    grouped={c:0 for c in centroids}
    for member,center in members.items():
        assert member in inputs and center in centroids
        grouped[center]+=inputs[member]
    assert grouped==centroids
    actual,found=set(truth),set(centroids)
    assert len(actual-found)==row['fn'] and len(found-actual)==row['fp']
    correct=wrong=unassigned=0
    origins={}
    with (source/'labels.csv').open() as stream:
        for label in csv.DictReader(stream):
            sequence=label['member']; count=int(label['read_count'])
            origins[sequence]=origins.get(sequence,0)+count
            center=members.get(sequence)
            if center is None:unassigned+=count
            elif center==label['true_barcode']:correct+=count
            else:wrong+=count
    assert origins==inputs
    assert correct+wrong+unassigned==total
    assert (correct,wrong,unassigned)==(row['correct_reads'],row['incorrect_reads'],row['unassigned_reads'])
    assert sum(centroids.values())+unassigned==total
    tv=(sum(abs(centroids.get(b,0)-truth.get(b,0)) for b in actual|found)+unassigned)/(2*total)
    close(tv,row['abundance_total_variation'])


def main():
    check_freeze()
    validation=load(HERE/'execution_validation.json')
    rows=load(HERE/'results.json')
    assert len(rows)==validation['cells']==726
    addendum=load(HERE/'reporting_addendum.json')
    for name,digest in addendum['files'].items():assert sha(HERE/name)==digest
    for r in rows:
        if r['status']!='complete':continue
        assert r['tp']+r['fn']==r['true_barcodes']
        assert r['tp']+r['fp']==r['inferred_barcodes']
        close(r['f1_percent'],200*r['tp']/(r['true_barcodes']+r['inferred_barcodes']))
        assert r['correct_reads']+r['incorrect_reads']+r['unassigned_reads']==r['input_reads']
        close(r['read_assignment_accuracy_percent'],100*r['correct_reads']/r['input_reads'])
    metric_names=['fn','fp','f1_percent','positive_truth_f1_percent','incorrect_reads','unassigned_reads',
                  'read_assignment_accuracy_percent','abundance_total_variation','input_reads','zero_read_truth']
    for filename in ['summary.json','successful_only_summary.json']:
        for r in load(HERE/filename):
            cells=[c for c in rows if c['condition']==r['condition'] and c['method']==r['method']]
            good=[c for c in cells if c['status']=='complete']
            assert len(good)==r['n_successful'] and len(cells)-len(good)==r['n_failed']
            if filename=='summary.json' and not r['complete_population']:
                assert all(r[m] is None for m in metric_names+['workflow_seconds'])
                continue
            if not good:continue
            for metric in metric_names:
                values=[c[metric] for c in good]
                close(r[metric],statistics.mean(values))
                if len(values)>1:close(r[metric+'_sd'],statistics.stdev(values))
            close(r['workflow_seconds'],statistics.median(c['workflow_seconds'] for c in good))
    contrasts=load(HERE/'accuracy_contrasts.json')
    dest=HERE/'generated/independent_validation'
    dest.mkdir(exist_ok=True)
    with (dest/'differences.csv').open('w') as stream:
        writer=csv.writer(stream);writer.writerow(['contrast','delta'])
        for i,c in enumerate(contrasts,1):
            if c['status']!='complete':continue
            left={r['seed']:r for r in rows if r['condition']==c['condition'] and r['method']=='lv'}
            right={r['seed']:r for r in rows if r['condition']==c['condition'] and r['method']==c['competitor']}
            assert left.keys()==right.keys() and len(left)==60
            for seed in sorted(left):writer.writerow([i,left[seed]['f1_percent']-right[seed]['f1_percent']])
    r_code='''args <- commandArgs(TRUE)
d <- read.csv(args[1])
out <- do.call(rbind, lapply(1:8, function(i) {
 x <- d$delta[d$contrast == i]; n <- length(x)
 if (n == 0) return(data.frame(contrast=i,n=0,mean=NA,sd=NA,lower=NA,p=NA))
 stopifnot(n == 60); m <- mean(x); s <- sd(x); se <- s/sqrt(n)
 p <- if (se == 0) as.numeric(m <= 0) else pt(m/se,n-1,lower.tail=FALSE)
 data.frame(contrast=i,n=n,mean=m,sd=s,lower=m-qt(1-.05/8,n-1)*se,p=p)
}))
out$holm <- p.adjust(ifelse(is.na(out$p),1,out$p),method="holm")
write.csv(out,args[2],row.names=FALSE,na="")
'''
    subprocess.run(['Rscript','-e',r_code,str(dest/'differences.csv'),str(dest/'r_inference.csv')],check=True)
    with (dest/'r_inference.csv').open() as stream:
        for observed in csv.DictReader(stream):
            expected=contrasts[int(observed['contrast'])-1]
            if expected['status']!='complete':
                assert int(observed['n'])==0
                continue
            for left,right in [('mean','mean_difference'),('sd','sd_difference'),('lower','simultaneous_one_sided_lower'),('p','one_sided_p'),('holm','holm_adjusted_p')]:
                close(float(observed[left]),expected[right])
    selected=load(HERE/'timing_results.json')
    original_lookup={key(r):r for r in rows}
    repairs={r['key']:r for r in load(HERE/'timing_repair_results.json')}
    assert len(selected)==726 and len(repairs)==28
    for row in selected:
        original=original_lookup[key(row)]
        for metric in ['fn','fp','f1_percent','correct_reads','incorrect_reads','unassigned_reads']:
            assert row.get(metric)==original.get(metric)
        if key(row) in repairs:
            repair=repairs[key(row)]
            assert repair['status']=='complete', 'Unsuccessful timing repair remains unavailable'
            close(row['selected_workflow_seconds'],repair['workflow_seconds'])
        elif original['status']=='complete':
            close(row['selected_workflow_seconds'],original['workflow_seconds'])
        else:assert row['selected_workflow_seconds'] is None
    for summary_row in load(HERE/'publication_summary.json'):
        timing_cells=[r for r in selected if r['condition']==summary_row['condition'] and r['method']==summary_row['method']]
        times=[r['selected_workflow_seconds'] for r in timing_cells if r['selected_workflow_seconds'] is not None]
        if len(times)==summary_row['planned_n']:
            close(summary_row['workflow_seconds'],statistics.median(times))
        else:assert summary_row['workflow_seconds'] is None
    timing_contrasts=load(HERE/'timing_contrasts.json')
    with (dest/'log_times.csv').open('w') as stream:
        writer=csv.writer(stream);writer.writerow(['contrast','delta'])
        for i,c in enumerate(timing_contrasts,1):
            if c['status']!='complete':continue
            left={r['seed']:r for r in selected if r['condition']==c['condition'] and r['method']=='lv'}
            right={r['seed']:r for r in selected if r['condition']==c['condition'] and r['method']==c['competitor']}
            for seed in sorted(left):writer.writerow([i,math.log(left[seed]['selected_workflow_seconds']/right[seed]['selected_workflow_seconds'])])
    timing_r='''args <- commandArgs(TRUE); d <- read.csv(args[1])
out <- do.call(rbind,lapply(1:8,function(i) {
 x <- d$delta[d$contrast == i]
 if (length(x)==0) return(data.frame(contrast=i,ratio=NA,lower=NA,upper=NA))
 stopifnot(length(x)==60); m <- mean(x); h <- qt(.975,59)*sd(x)/sqrt(60)
 data.frame(contrast=i,ratio=exp(m),lower=exp(m-h),upper=exp(m+h))
}))
write.csv(out,args[2],row.names=FALSE,na="")
'''
    subprocess.run(['Rscript','-e',timing_r,str(dest/'log_times.csv'),str(dest/'r_timing.csv')],check=True)
    with (dest/'r_timing.csv').open() as stream:
        for observed in csv.DictReader(stream):
            expected=timing_contrasts[int(observed['contrast'])-1]
            if expected['status']!='complete':continue
            for left,right in [('ratio','geometric_time_ratio'),('lower','ci95_lower'),('upper','ci95_upper')]:
                close(float(observed[left]),expected[right])
    campaign_start=min(datetime.fromisoformat(r['started_at']) for r in rows)
    campaign_end=max(datetime.fromisoformat(r['finished_at']) for r in repairs.values())
    power_log=subprocess.check_output(['pmset','-g','log'],text=True)
    sleeps=[]
    for line in power_log.splitlines():
        if 'Entering Sleep state' not in line:continue
        when=datetime.strptime(line[:25],'%Y-%m-%d %H:%M:%S %z')
        if campaign_start<=when<=campaign_end:
            assert START<=when<=END, 'Additional unregistered sleep interval: '+line
            sleeps.append(line)
    expected_sleeps=[s for s in load(HERE/'timing_repair_protocol.json')['power_log_excerpt'] if 'Entering Sleep state' in s]
    assert len(sleeps)==len(expected_sleeps)==4, 'Power log coverage must retain every known sleep event'
    save(HERE/'timing_environment_validation.json',dict(no_additional_sleep_intervals=True,
        registered_interval_sleep_events=sleeps,first_call=campaign_start.isoformat(),last_repair=campaign_end.isoformat()))
    datasets=load(HERE/'datasets.json')
    first=load(HERE/'final_protocol.json')['final_seeds'][0]
    keys={(first,'random_mixed'),(first,'anchored_mixed')}
    largest=max(datasets,key=lambda d:d['event_counts'].get('boundary_reads',0))
    keys.add((largest['seed'],largest['condition']))
    selected=[r for r in rows if r['status']=='complete' and (r['condition']=='milos' or (r['seed'],r['condition']) in keys)]
    for row in selected:
        key='milos' if row['condition']=='milos' else f"{row['seed']}/{row['condition']}"
        source=REFERENCE/'generated/datasets/milos' if key=='milos' else HERE/'generated/datasets'/key
        check_mapping(source,HERE/'generated/results'/key/row['method'],row)
    save(HERE/'independent_validation.json',dict(registered_rows=len(rows),successful_row_arithmetic_checked=validation['successful_cells'],
        every_summary_mean_sd_and_median_checked=True,available_paired_t_and_holm_checked_with_base_R=True,
        sampled_mapping_cells=len(selected),sample_selection='All successful methods on first registered seed in both designs, highest boundary-read input, and fixed Milo',
        sampled_mapping_counts_and_origins_reconcile=True,bootstrap_scope='Implementation and frozen-seed checks; no independent bootstrap implementation',
        reporting_addendum_hashes_match=True,original_frozen_reference_preserved=True,
        selected_timing_medians_and_log_time_intervals_checked=True,
        timing_intervals_independently_checked_with_base_R=True,
        sleep_affected_observations_repaired_without_selecting_minimum_time=True))
    print('INDEPENDENT VALIDATION PASSED',len(selected),'mapping samples')


if __name__=='__main__':
    main()
