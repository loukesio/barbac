"""Serial, label-blind comparison on the existing 100k reference simulation.

Large inputs/outputs stay in an ignored work directory. Each completed method
is checkpointed; --resume evaluates existing outputs without rerunning tools.
Dependencies: pandas, numpy, rapidfuzz; R barbac/readr; external method binaries.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np
import pandas as pd
from rapidfuzz import process
from rapidfuzz.distance import Levenshtein

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def profile(source):
    inputs = pd.read_csv(source / 'barbac_benchmark_input.csv')
    labels = pd.read_csv(source / 'simulated_reads.csv')
    truth = pd.read_csv(source / 'true_counts.csv').rename(columns={'BC':'barcode', 'True Count':'true_count'})
    assert list(inputs.columns) == ['barcode', 'counts']
    assert list(labels.columns) == ['BC', 'Count', 'true_BC']
    for frame, key, count in [(inputs,'barcode','counts'),(labels,'BC','Count'),(truth,'barcode','true_count')]:
        assert not frame.isna().any().any()
        assert frame[key].is_unique
        assert frame[key].str.fullmatch('[ACGT]+').all()
        assert (frame[count] >= (0 if frame is truth else 1)).all()
        assert (frame[count] % 1 == 0).all()
    assert labels.true_BC.isin(truth.barcode).all()
    observed = labels.set_index('BC').Count.sort_index()
    pd.testing.assert_series_equal(inputs.set_index('barcode').counts.sort_index(), observed, check_names=False)
    realized = labels.groupby('true_BC').Count.sum().reindex(truth.barcode, fill_value=0)
    # Preserve source inconsistencies rather than silently rewriting labels.
    delta = realized.to_numpy() - truth.true_count.to_numpy()
    assert int(inputs.counts.sum()) == int(truth.true_count.sum())
    mismatches = truth.loc[delta != 0].copy()
    mismatches['label_count'] = realized.to_numpy()[delta != 0]
    lengths = inputs.barcode.str.len()
    info = dict(input_sequences=len(inputs), true_barcodes=len(truth), input_reads=int(inputs.counts.sum()),
                absent_truth=int((~truth.barcode.isin(inputs.barcode)).sum()),
                zero_read_truth=int((truth.true_count == 0).sum()),
                parent_count_mismatches=mismatches.to_dict(orient='records'),
                parent_count_absolute_difference=int(np.abs(delta).sum()),
                length_distribution={str(k):int(v) for k,v in lengths.value_counts().items()},
                off_modal_reads=int(inputs.loc[lengths != lengths.mode()[0], 'counts'].sum()),
                checks='PASS: unique nonnull ACGT keys; positive integer observed counts/nonnegative truth counts; exact input and source counts; valid parent keys; grand totals reconcile. Per-parent differences reported separately.',
                source_sha256={name:sha(source/name) for name in ['barbac_benchmark_input.csv','simulated_reads.csv','true_counts.csv']})
    return inputs, labels, truth, info


def evaluate(centroids, members, labels, truth):
    assert list(centroids.columns) == ['central_barcode','sum_counts']
    assert not centroids.isna().any().any()
    assert (centroids.sum_counts > 0).all()
    assert members.member.is_unique and not members.isna().any().any()
    assert members.member.isin(labels.BC).all()
    assert members.central_barcode.isin(centroids.central_barcode).all()
    joined = labels.merge(members[['member','central_barcode']], left_on='BC', right_on='member', how='left', validate='one_to_one')
    expected_counts = joined.dropna(subset=['central_barcode']).groupby('central_barcode').Count.sum().sort_index()
    actual_counts = centroids.groupby('central_barcode').sum_counts.sum().sort_index()
    pd.testing.assert_series_equal(actual_counts, expected_counts, check_names=False, check_dtype=False)
    truth_set, found = set(truth.barcode), set(centroids.central_barcode)
    absent = truth_set - set(labels.BC)
    missing, extra = truth_set-found, found-truth_set
    tp, fn, fp = len(truth_set & found), len(missing), len(extra)
    positive_truth = set(truth.loc[truth.true_count > 0, 'barcode'])
    positive_tp = len(positive_truth & found)
    positive_fp = len(found - positive_truth)
    positive_fn = len(positive_truth - found)
    total = int(labels.Count.sum())
    correct = int(joined.loc[joined.central_barcode == joined.true_BC, 'Count'].sum())
    unassigned = int(joined.loc[joined.central_barcode.isna(), 'Count'].sum())
    realized = labels.groupby('true_BC').Count.sum().reindex(truth.barcode,fill_value=0).to_numpy()
    inconsistent = set(truth.loc[realized != truth.true_count.to_numpy(), 'barcode'])
    consistent = ~joined.true_BC.isin(inconsistent)
    consistent_total = int(joined.loc[consistent, 'Count'].sum())
    consistent_correct = int(joined.loc[consistent & (joined.central_barcode == joined.true_BC), 'Count'].sum())
    ws = 0
    extras = sorted(extra)
    for start in range(0,len(extras),100):
        distances = process.cdist(extras[start:start+100], sorted(truth_set), scorer=Levenshtein.distance, score_cutoff=3, dtype=np.uint8)
        ws += int((distances.min(axis=1)<=3).sum())
    return dict(centroids=len(centroids), unique_centroids=len(found), tp=tp, fn=fn, fp=fp,
                observed_fn=len(missing-absent), absent_truth_recovered=len(found & absent), ws=ws,
                precision=tp/(tp+fp), recall=tp/len(truth_set), f1=2*tp/(2*tp+fn+fp),
                positive_truth_fn=positive_fn, positive_truth_fp=positive_fp,
                positive_truth_f1=2*positive_tp/(2*positive_tp+positive_fn+positive_fp),
                correct_reads=correct, misassigned_reads=total-correct-unassigned, unassigned_reads=unassigned,
                read_assignment_accuracy=correct/total, output_reads=int(centroids.sum_counts.sum()),
                consistent_parent_reads=consistent_total,
                consistent_parent_accuracy=consistent_correct/consistent_total,
                input_reads=total, mapping_counts_reconcile=True)


def execute(command, work, log, extra_env=None):
    start=time.perf_counter()
    with log.open('w') as stream:
        subprocess.run([str(x) for x in command], cwd=work, stdout=stream, stderr=subprocess.STDOUT, check=True,
                       env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1', **(extra_env or {})})
    elapsed = time.perf_counter()-start
    log.with_suffix('.process.json').write_text(json.dumps({'command':[str(x) for x in command], 'process_seconds':elapsed},indent=2)+'\n')
    return elapsed


def run_method(name, args, work, inputs):
    dest=work/name
    dest.mkdir(exist_ok=True)
    checkpoint=dest/'run.json'
    if args.resume and checkpoint.exists():
        saved = json.loads(checkpoint.read_text())
        for filename, expected in saved['output_sha256'].items():
            assert sha(dest/filename) == expected, (name, filename)
        return saved
    start=time.perf_counter()
    prep=0.0
    core=None
    build=None
    if name.startswith(('barbac_','previous_')):
        old=name.startswith('previous_')
        method='hamming' if 'hamming' in name else 'lv'
        ordering='support' if name.endswith('support') else 'sequence'
        lib=args.baseline_library if old else args.library
        command=[args.rscript,HERE/'run_barbac.R',lib,work/'input.csv',dest,method,ordering]
        wall=execute(command,dest,dest/'run.log')
        for line in (dest/'run.log').read_text().splitlines():
            if line.startswith('BARBAC_ALGO_SECONDS='): core=float(line.split('=')[1])
            if line.startswith('BARBAC_BUILD_ID='): build=line.split('=')[1]
        assert core is not None and build is not None
    elif name=='shepherd':
        command=[sys.executable,args.tools/'Shepherd/shepherd_t0.py','-f',work/'input.tsv','-l','20','-eps','3']
        wall=execute(command,dest,dest/'run.log')
        # Shepherd derives output paths from its absolute input prefix.
        c=pd.read_csv(work/'input_pb_freq.csv');c.columns=['central_barcode','sum_counts']
        m=pd.read_csv(work/'input_seq_clust.csv')
        own=m[m.sequence.isin(c.central_barcode)]
        assert own.cluster.is_unique and len(own)==len(c)
        mapping=own.set_index('cluster').sequence
        m=pd.DataFrame({'member':m.sequence,'central_barcode':m.cluster.map(mapping)})
        # Shepherd emits a -1 cluster for discarded reads on some datasets.
        m=m.dropna(subset=['central_barcode'])
        c.to_csv(dest/'centroids.csv',index=False);m.to_csv(dest/'members.csv',index=False)
    elif name.startswith('starcode_'):
        command=[args.tools/'starcode/starcode','-d','3','-t','1','--print-clusters','-i',work/'input.tsv','-o',dest/'clusters.tsv']
        if name=='starcode_sphere': command+=['-s']
        wall=execute(command,dest,dest/'run.log')
        raw=pd.read_csv(dest/'clusters.tsv',sep='\t',header=None,names=['central_barcode','sum_counts','members'])
        raw[['central_barcode','sum_counts']].to_csv(dest/'centroids.csv',index=False)
        m=raw.assign(member=raw.members.str.split(',')).explode('member')[['member','central_barcode']]
        m.to_csv(dest/'members.csv',index=False)
    elif name=='bartender':
        t0=time.perf_counter()
        with (dest/'expanded.csv').open('w') as out:
            read_id=0
            for seq,cnt in inputs.itertuples(index=False,name=None):
                out.write(''.join(f'{seq},{i}\n' for i in range(read_id+1,read_id+int(cnt)+1)))
                read_id+=int(cnt)
        prep=time.perf_counter()-t0
        command=[args.tools/'bartender-1.1/bartender_single_com','-f',dest/'expanded.csv','-o',dest/'bartender','-d','3','-t','1']
        wall=execute(command,dest,dest/'run.log', {'PATH':str(args.tools/'bartender-1.1')+os.pathsep+os.environ.get('PATH','')})
        c=pd.read_csv(dest/'bartender_cluster.csv')
        raw=pd.read_csv(dest/'bartender_barcode.csv')
        m=pd.DataFrame({'member':raw['Unique.reads'],'central_barcode':raw['Cluster.ID'].map(c.set_index('Cluster.ID').Center)})
        c.rename(columns={'Center':'central_barcode','time_point_1':'sum_counts'})[['central_barcode','sum_counts']].to_csv(dest/'centroids.csv',index=False)
        m.to_csv(dest/'members.csv',index=False)
        (dest/'expanded.csv').unlink()
    else: raise ValueError(name)
    result=dict(method=name,process_seconds=wall,core_seconds=core,input_preparation_seconds=prep,
                pipeline_seconds=time.perf_counter()-start,build_id=build,command=[str(x) for x in command],
                output_sha256={f:sha(dest/f) for f in ['centroids.csv','members.csv']})
    checkpoint.write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--tools',type=Path,required=True)
    p.add_argument('--library',type=Path,required=True)
    p.add_argument('--baseline-library',type=Path)
    p.add_argument('--rscript',default='/usr/local/bin/Rscript')
    p.add_argument('--work',type=Path,default=ROOT/'benchmark/four_condition_comparison/generated/reference_2026-09-08')
    p.add_argument('--report',type=Path,default=HERE)
    p.add_argument('--resume',action='store_true')
    p.add_argument('--methods',nargs='+',default=['barbac_hamming_sequence','barbac_hamming_support','barbac_lv_sequence','barbac_lv_support','previous_hamming_sequence','previous_lv_sequence','shepherd','starcode_sphere','starcode_mp','bartender'])
    args=p.parse_args()
    if any(name.startswith('previous_') for name in args.methods) and args.baseline_library is None:
        p.error('--baseline-library is required when selecting previous methods')
    for key in ['source','tools','library','work','report']:
        setattr(args,key,getattr(args,key).resolve())
    args.work.mkdir(parents=True,exist_ok=True);args.report.mkdir(parents=True,exist_ok=True)
    inputs,labels,truth,quality=profile(args.source)
    (args.report/'data_quality.json').write_text(json.dumps(quality,indent=2)+'\n')
    # No true_BC/truth counts enter the sorting or any tool input.
    inputs=inputs.sort_values(['counts','barcode'],ascending=[False,True]).reset_index(drop=True)
    if not args.resume or not (args.work/'input.csv').exists():
        inputs.to_csv(args.work/'input.csv',index=False)
        inputs.to_csv(args.work/'input.tsv',sep='\t',index=False,header=False)
    else:
        pd.testing.assert_frame_equal(pd.read_csv(args.work/'input.csv'),inputs)
    provenance=dict(platform=platform.platform(),processor=platform.processor(),python=sys.version,
                    python_packages={name:importlib.metadata.version(name) for name in ['pandas','numpy','scipy','rapidfuzz']},
                    r_version=subprocess.check_output([args.rscript,'--version'],text=True).strip(),
                    git_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                    input_sha256={f:sha(args.work/f) for f in ['input.csv','input.tsv']},
                    source_sha256={str(x.relative_to(ROOT)):sha(x) for x in [HERE/'run_comparison.py',HERE/'run_barbac.R',ROOT/'src/clustering.cpp',ROOT/'R/11_super_cluster2.R']},
                    libraries={str(lib):{str(f.relative_to(lib)):sha(f) for f in (lib/'barbac').rglob('*') if f.is_file()} for lib in [args.library,args.baseline_library] if lib},
                    tools={name:subprocess.check_output(['git','rev-parse','HEAD'],cwd=args.tools/name,text=True).strip() for name in ['Shepherd','starcode','bartender-1.1']},
                    tool_file_sha256={name:sha(args.tools/name) for name in ['Shepherd/shepherd_t0.py','starcode/starcode','bartender-1.1/bartender_single_com','bartender-1.1/bartender_single']},
                    protocol='One serial run per configuration; fixed abundance then sequence ordering; distance 3; barbac ratio20/error0.005; single-thread environment. Core includes CSV input; process includes native exports/startup; pipeline also includes format conversion. Common input staging and metric evaluation excluded. No labels used for clustering.')
    (args.report/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    results=[]
    for name in args.methods:
        print(f'START {name}',flush=True)
        run=run_method(name,args,args.work,inputs)
        dest=args.work/name
        scores=evaluate(pd.read_csv(dest/'centroids.csv'),pd.read_csv(dest/'members.csv'),labels,truth)
        result={**{k:v for k,v in run.items() if k not in ['command','output_sha256']},**scores}
        results.append(result)
        pd.DataFrame(results).to_csv(args.report/'results.csv',index=False)
        (args.report/f'{name}.json').write_text(json.dumps({**run,**scores},indent=2)+'\n')
        print(json.dumps(result),flush=True)

if __name__=='__main__': main()
