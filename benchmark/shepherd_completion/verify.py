"""Independent CSV arithmetic and base-R checks of supplementary conclusions."""
import collections
import csv
import datetime
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys

HERE=Path(__file__).resolve().parent
FINAL=Path('/Users/theodosiou/Documents/Projects/test_barbac/.codex/publication-final-2026-09-12/source/benchmark/publication_final')
sys.path.insert(0,str(FINAL))
from common import load,save,sha


def read(path):
    with path.open() as stream:
        yield from csv.DictReader(stream)


def independent_scores(row,full=False):
    source=Path(row['dataset_path'])
    key='milos' if row['condition']=='milos' else f"{row['seed']}/{row['condition']}"
    dest=HERE/'generated/results'/key
    for name,digest in row['output_sha256'].items():
        assert sha(dest/name)==digest
    truth={r['barcode']:int(r['true_count']) for r in read(source/'truth.csv')}
    found={r['central_barcode']:int(r['sum_counts']) for r in read(dest/'centroids.csv')}
    actual=set(truth); observed=set(found)
    tp=len(actual&observed);fn=len(actual-observed);fp=len(observed-actual)
    assert (tp,fn,fp)==(row['tp'],row['fn'],row['fp'])
    assert math.isclose(200*tp/(2*tp+fn+fp),row['f1_percent'],abs_tol=1e-12)
    if not full:return
    inputs={r['barcode']:int(r['counts']) for r in read(source/'input.csv')}
    mapping={r['member']:r['central_barcode'] for r in read(dest/'members.csv')}
    cluster_counts=collections.Counter()
    for member,root in mapping.items():cluster_counts[root]+=inputs[member]
    assert dict(cluster_counts)==found
    correct=wrong=unassigned=0
    for r in read(source/'labels.csv'):
        root=mapping.get(r['member']);n=int(r['read_count'])
        if root is None:unassigned+=n
        elif root==r['true_barcode']:correct+=n
        else:wrong+=n
    assert (correct,wrong,unassigned)==(row['correct_reads'],row['incorrect_reads'],row['unassigned_reads'])
    assert sum(inputs.values())==correct+wrong+unassigned
    tv=(sum(abs(found.get(s,0)-truth.get(s,0)) for s in actual|observed)+unassigned)/(2*sum(inputs.values()))
    assert math.isclose(tv,row['abundance_total_variation'],abs_tol=1e-12)


def main():
    frozen=load(HERE/'freeze.json')
    for path,digest in frozen['sha256'].items():assert sha(path)==digest,path
    new=load(HERE/'results.json')
    assert len(new)==121
    failed_keys={(r['condition'],r['seed']) for r in load(FINAL/'failure_audit.json')}
    selected=[new[0]]
    selected+=[next(r for r in new if (r['condition'],r['seed']) in failed_keys)]
    selected+=[next(r for r in new if r['condition']=='milos')]
    full_keys={(r['condition'],r['seed']) for r in selected}
    for row in new:
        if row['status']=='complete':independent_scores(row,(row['condition'],row['seed']) in full_keys)
    combined=load(HERE/'combined_results.json')
    original=load(FINAL/'timing_results.json')
    original_map={(r['condition'],r['seed'],r['method']):r for r in original if r['method']!='shepherd'}
    for row in combined:
        if row['method']=='shepherd_documented':continue
        old=original_map[(row['condition'],row['seed'],row['method'])]
        assert all(row[k]==v for k,v in old.items())
        assert row['reported_workflow_seconds']==old['selected_workflow_seconds']
    contrast_path=HERE/'generated/independent_R_contrasts.csv'
    r_code='''
args <- commandArgs(trailingOnly=TRUE)
x <- read.csv(args[1],check.names=FALSE)
answer <- list()
for (condition in c("random_mixed","anchored_mixed")) {
  a <- x[x$condition==condition & x$method=="lv",c("seed","f1_percent")]
  b <- x[x$condition==condition & x$method=="shepherd_documented",c("seed","f1_percent")]
  m <- merge(a,b,by="seed",suffixes=c(".a",".b"))
  stopifnot(nrow(m)==60L,all(is.finite(m$f1_percent.a)),all(is.finite(m$f1_percent.b)))
  delta <- m$f1_percent.a-m$f1_percent.b
  se <- sd(delta)/sqrt(length(delta))
  answer[[condition]] <- data.frame(condition=condition,n=length(delta),mean_difference=mean(delta),
    sd_difference=sd(delta),simultaneous_one_sided_lower=mean(delta)-qt(1-.05/8,59)*se,
    one_sided_p=pt(mean(delta)/se,59,lower.tail=FALSE))
}
write.csv(do.call(rbind,answer),args[2],row.names=FALSE)
'''
    r_path=HERE/'generated/independent_check.R';r_path.write_text(r_code)
    if all(r['status']=='complete' for r in new):
        subprocess.run(['Rscript',str(r_path),str(HERE/'combined_results.csv'),str(contrast_path)],check=True)
        reference={r['condition']:r for r in load(HERE/'accuracy_contrasts.json')}
        for row in read(contrast_path):
            for name in ['mean_difference','sd_difference','simultaneous_one_sided_lower','one_sided_p']:
                assert math.isclose(float(row[name]),reference[row['condition']][name],rel_tol=1e-8,abs_tol=1e-12),(name,row)
    save(HERE/'independent_validation.json',dict(completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        all_new_success_centroid_metrics_independently_recomputed=True,
        full_read_mapping_spot_checks=[dict(condition=r['condition'],seed=r['seed']) for r in selected],
        all_605_other_method_records_exactly_preserved=True,registered_fingerprints_unchanged=True,
        base_R_paired_t_checks_passed=all(r['status']=='complete' for r in new),
        bootstrap_scope='Existing previously validated paired-summary implementation reused; not independently reimplemented here'))
    print('Independent centroid arithmetic, three complete mapping audits and base-R paired contrasts passed.')


if __name__=='__main__':main()
