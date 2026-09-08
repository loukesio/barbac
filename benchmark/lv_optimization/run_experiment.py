"""Serial ablation: v12 LV, v13 exact pruning, v13 Poisson indel exception.

Reuses hashed reference and four-condition inputs. Large outputs are ignored.
Run with the Python environment containing pandas/numpy/rapidfuzz/scipy.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(ROOT/'benchmark/reference_comparison'))
from run_comparison import evaluate, profile, sha


def run(library, source, output, tie, model, rate=0.005):
    output.mkdir(parents=True, exist_ok=True)
    command=['/usr/local/bin/Rscript',str(HERE/'run_barbac.R'),str(library),str(source),str(output),'lv',tie]
    if model is not None: command += [model, str(rate)]
    start=time.perf_counter()
    with (output/'run.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True,
                       env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1'})
    wall=time.perf_counter()-start
    markers=dict(line.split('=',1) for line in (output/'run.log').read_text().splitlines()
                 if line.startswith(('BARBAC_ALGO_SECONDS=','BARBAC_BUILD_ID=')))
    return {'core_seconds':float(markers['BARBAC_ALGO_SECONDS']),'process_seconds':wall,
            'build_id':markers['BARBAC_BUILD_ID'],'command':command,
            'output_sha256':{f:sha(output/f) for f in ['centroids.csv','members.csv']}}


def validate_counts(input_path, output):
    data=pd.read_csv(input_path).set_index('barcode').counts
    members=pd.read_csv(output/'members.csv')
    centroids=pd.read_csv(output/'centroids.csv')
    assert members.member.is_unique and centroids.central_barcode.is_unique
    actual=members.set_index('member').member_count.sort_index()
    pd.testing.assert_series_equal(data.sort_index(),actual,check_names=False)
    pd.testing.assert_series_equal(members.groupby('central_barcode').member_count.sum().sort_index(),
                                  centroids.set_index('central_barcode').sum_counts.sort_index(),check_names=False)
    return centroids


def compare_centroids(output, previous):
    a=pd.read_csv(output/'centroids.csv').sort_values('central_barcode').reset_index(drop=True)
    b=pd.read_csv(previous).sort_values('central_barcode').reset_index(drop=True)
    pd.testing.assert_frame_equal(a,b)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--library',type=Path,default=Path('/private/tmp/barbac-lv-v13'))
    parser.add_argument('--baseline-library',type=Path,default=Path('/private/tmp/barbac-exact-final'))
    parser.add_argument('--source',type=Path,default=Path('/Users/theodosiou/Documents/Projects/Barcodes/barbac-benchmark'))
    parser.add_argument('--work',type=Path,default=ROOT/'benchmark/four_condition_comparison/generated/lv-v13')
    parser.add_argument('--phase',choices=['reference','four','all'],default='all')
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    args.work.mkdir(parents=True,exist_ok=True)
    sources=['src/clustering.cpp','R/11_super_cluster2.R','R/RcppExports.R','src/RcppExports.cpp','benchmark/lv_optimization/run_barbac.R','benchmark/lv_optimization/run_experiment.py']
    manifest={'platform':platform.platform(),'source_sha256':{f:sha(ROOT/f) for f in sources},
              'libraries':{str(lib):{str(f.relative_to(lib)):sha(f) for f in (lib/'barbac').rglob('*') if f.is_file()} for lib in [args.library,args.baseline_library]},
              'configuration':{'distance':3,'merge_ratio':20,'error_rate':0.005,'poisson_upper_tail_threshold':0.01,'phase':args.phase},
              'timing':'Serial fresh R processes, core includes input/sorting/clustering, process includes startup and both exports. Scoring excluded.'}
    manifest_path=args.work/'manifest.json'
    if args.resume and manifest_path.exists():
        saved=json.loads(manifest_path.read_text())
        assert saved['source_sha256']==manifest['source_sha256'] and saved['libraries']==manifest['libraries']
    manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
    rows=[]
    def measure(dataset,seed,variant,tie,lib,model,source,truth,old_output=None,labels=None):
        output=args.work/dataset/str(seed)/variant/tie
        checkpoint=output/'result.json'
        if args.resume and checkpoint.exists():
            result=json.loads(checkpoint.read_text())
            for f,h in result['output_sha256'].items(): assert sha(output/f)==h
        else:
            print('START',dataset,seed,variant,tie,flush=True)
            result=run(lib,source,output,tie,model)
            centroids=validate_counts(source,output)
            if labels is not None:
                scores=evaluate(centroids,pd.read_csv(output/'members.csv'),labels,truth)
            else:
                true_set=set(truth.iloc[:,0]);found=set(centroids.central_barcode)
                tp=len(true_set&found);fn=len(true_set-found);fp=len(found-true_set)
                scores={'fn':fn,'fp':fp,'f1':2*tp/(2*tp+fn+fp),'output_reads':int(centroids.sum_counts.sum())}
            if old_output is not None:
                compare_centroids(output,old_output)
                result['centroid_equivalence']='passed'
                prior_members=old_output.parent/'members.csv'
                if prior_members.exists():
                    a=pd.read_csv(output/'members.csv').sort_values('member').reset_index(drop=True)
                    b=pd.read_csv(prior_members).sort_values('member').reset_index(drop=True)
                    pd.testing.assert_frame_equal(a,b)
                    result['membership_equivalence']='passed'
            result.update(scores,dataset=dataset,seed=seed,variant=variant,tie_break=tie,input_sha256=sha(source))
            checkpoint.write_text(json.dumps(result,indent=2)+'\n')
        rows.append({k:v for k,v in result.items() if k not in ['command','output_sha256']})
        pd.DataFrame(rows).to_csv(HERE/f'results_{args.phase}.csv',index=False)
        print(json.dumps(rows[-1]),flush=True)
    if args.phase in ['all','reference']:
        _,labels,truth,quality=profile(args.source)
        source=ROOT/'benchmark/four_condition_comparison/generated/reference_2026-09-08/input.csv'
        original=json.loads((ROOT/'benchmark/reference_comparison/provenance.json').read_text())
        assert sha(source)==original['input_sha256']['input.csv']
        for tie in ['support','sequence']:
            old=source.parent/f'barbac_lv_{tie}'/'centroids.csv'
            measure('reference',0,'pruned',tie,args.library,'none',source,truth,old,labels)
            measure('reference',0,'poisson',tie,args.library,'poisson',source,truth,labels=labels)
    if args.phase in ['all','four']:
        for seed in [42,43,44]:
            for condition in ['random_substitutions','random_low_indels','anchored_substitutions','anchored_low_indels']:
                directory=ROOT/'benchmark/four_condition_comparison/generated/revision-accuracy'/str(seed)/condition
                saved=json.loads((directory/'simulation.json').read_text())
                for f,h in saved['sha256'].items(): assert sha(directory/f)==h
                source=directory/'input.csv';truth=pd.read_csv(directory/'true_counts.csv')
                for variant,lib,model in [('v12',args.baseline_library,None),('pruned',args.library,'none'),('poisson',args.library,'poisson')]:
                    old=directory/'candidate_support_lv_0.csv' if variant in ['v12','pruned'] else None
                    measure(condition,seed,variant,'support',lib,model,source,truth,old)
    (HERE/f'manifest_{args.phase}.json').write_text(json.dumps(manifest,indent=2)+'\n')

if __name__=='__main__': main()
