"""Fresh serial seed-42 peer runs, including membership exports and validation.

Run only after run_comparison.py finishes, to avoid concurrent clustering.
"""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import pandas as pd
from run_comparison import HERE, ROOT, BASE, CONDITIONS, sha, score


def main():
    work = BASE/'generated/latest-four-v13/peers'
    assert json.loads((HERE/'manifest.json').read_text())['status']=='complete'
    tools = Path('/Users/theodosiou/Documents/Projects/Barcodes/barbac-benchmark/tools')
    prior = json.loads((BASE/'versions_2026-09-08.json').read_text())
    versions = {name:subprocess.check_output(['git','rev-parse','HEAD'],cwd=tools/name,text=True).strip()
                for name in ['Shepherd','starcode']}
    assert versions['Shepherd']==prior['shepherd_revision']
    assert versions['starcode']==prior['starcode_revision']
    env = {**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1',
           'MKL_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1','PYTHONHASHSEED':'0'}
    rows=[]
    receipts={}
    for condition in CONDITIONS:
        source=BASE/'generated/revision-accuracy/42'/condition
        inputs,truth=pd.read_csv(source/'input.csv'),pd.read_csv(source/'true_counts.csv')
        cfg=json.loads((source/'simulation.json').read_text())
        for file,expected in cfg['sha256'].items(): assert sha(source/file)==expected
        for method in ['shepherd','starcode']:
            dest=work/condition/method
            dest.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(source/'shepherd_input.txt',dest/'input.tsv')
            assert sha(dest/'input.tsv')==cfg['sha256']['shepherd_input.txt']
            checkpoint=dest/'result.json'
            if checkpoint.exists():
                row=json.loads(checkpoint.read_text())
                for file,expected in row['hashes'].items(): assert sha(dest/file)==expected
            else:
                print('START peer',condition,method,flush=True)
                if method=='shepherd':
                    command=[sys.executable,str(tools/'Shepherd/shepherd_t0.py'),'-f',str(dest/'input.tsv'),
                             '-l',str(cfg['configuration']['barcode_length']),'-eps','3']
                else:
                    command=[str(tools/'starcode/starcode'),'-d','3','-s','-t','1','--print-clusters',
                             '-i',str(dest/'input.tsv'),'-o',str(dest/'clusters.tsv')]
                start=time.perf_counter()
                with (dest/'run.log').open('w') as log:
                    subprocess.run(command,cwd=dest,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
                process=time.perf_counter()-start
                if method=='shepherd':
                    c=pd.read_csv(dest/'input_pb_freq.csv');c.columns=['central_barcode','sum_counts']
                    raw=pd.read_csv(dest/'input_seq_clust.csv')
                    roots=raw[raw.sequence.isin(c.central_barcode)]
                    assert roots.cluster.is_unique and len(roots)==len(c)
                    members=pd.DataFrame({'member':raw.sequence,'central_barcode':raw.cluster.map(roots.set_index('cluster').sequence)}).dropna()
                else:
                    raw=pd.read_csv(dest/'clusters.tsv',sep='\t',header=None,names=['central_barcode','sum_counts','members'])
                    c=raw[['central_barcode','sum_counts']]
                    members=raw.assign(member=raw.members.str.split(',')).explode('member')[['member','central_barcode']]
                members['member_count']=members.member.map(inputs.set_index('barcode').counts)
                c.to_csv(dest/'centroids.csv',index=False)
                members.to_csv(dest/'members.csv',index=False)
                workflow=time.perf_counter()-start
                row=dict(seed=42,condition=condition,method=method,distance=3,indel_model='not_applicable',
                         origin='fresh_peer',core_seconds=None,process_seconds=process,workflow_seconds=workflow,
                         command=command,hashes={file:sha(dest/file) for file in ['centroids.csv','members.csv']},
                         output_sha256=sha(dest/'centroids.csv'))
            c=pd.read_csv(dest/'centroids.csv');members=pd.read_csv(dest/'members.csv')
            assert members.member.is_unique and not members.isna().any().any()
            assert members.member.isin(inputs.barcode).all()
            assert members.central_barcode.isin(c.central_barcode).all()
            pd.testing.assert_series_equal(members.groupby('central_barcode').member_count.sum().sort_index(),
                                          c.set_index('central_barcode').sum_counts.sort_index(),check_names=False,check_dtype=False)
            metrics=score(c,inputs,truth)
            old=pd.read_csv(BASE/'generated/exact-search/42'/condition/f'fresh_{method}.csv')
            row['previous_centroid_count_equivalence']=c.sort_values('central_barcode').reset_index(drop=True).equals(old.sort_values('central_barcode').reset_index(drop=True))
            row.update(metrics)
            checkpoint.write_text(json.dumps(row,indent=2)+'\n')
            rows.append({k:v for k,v in row.items() if k not in ['command','hashes']})
            receipts[f'{condition}/{method}']=row
            pd.DataFrame(rows).to_csv(HERE/'fresh_peers.csv',index=False)
            print(f"DONE peer {condition}/{method}: FN={row['fn']} FP={row['fp']} workflow={row['workflow_seconds']:.2f}s equal_previous={row['previous_centroid_count_equivalence']}",flush=True)
    manifest=dict(status='complete',tool_revisions=versions,
                  source_sha256=sha(HERE/'time_peers.py'),
                  tool_sha256={name:sha(tools/name) for name in ['Shepherd/shepherd_t0.py','starcode/starcode']},
                  environment= {k:env[k] for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','PYTHONHASHSEED']},
                  protocol='Serial, no concurrent clustering. Seed 42 only. Distance 3, Starcode sphere with one thread. Workflow includes process and conversion to centroid/member exports; common input staging and metric scoring excluded.',
                  historical_timing_correction='The earlier three-seed peer accuracy campaign overlapped another accuracy campaign. Its timings are exploratory and must not be used for clean speed ratios. Use these fresh seed-42 workflow timings.',
                  runs=receipts)
    (HERE/'peer_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__': main()
