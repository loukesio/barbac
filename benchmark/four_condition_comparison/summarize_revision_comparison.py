"""Validate and compact the measured September 8 revision experiment.

The barbac and competitor campaigns ran separately against independently
regenerated inputs. Join only after verifying byte-identical simulated inputs.
"""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd

HERE=Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--barbac-dir',type=Path,required=True)
    p.add_argument('--competitor-dir',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,default=HERE)
    args=p.parse_args()
    a=pd.read_csv(args.barbac_dir/'results.csv')
    c=pd.read_csv(args.competitor_dir/'fresh_competitors.csv')
    manifest=json.loads((args.barbac_dir/'manifest.json').read_text())
    assert manifest['status']=='complete'
    cfg=manifest['configuration']
    expected_pairs={(s,condition['slug']) for s in cfg['seeds'] for condition in cfg['conditions']}
    assert set(zip(a.seed,a.condition))==expected_pairs
    assert set(zip(c.seed,c.condition))==expected_pairs
    assert len(a)==len(expected_pairs)*len(cfg['libraries'])*2*cfg['repeats']
    assert len(c)==len(expected_pairs)*2
    inputs={}
    competitor_outputs={}
    for seed,condition in sorted(expected_pairs):
        ad=args.barbac_dir/str(seed)/condition
        cd=args.competitor_dir/str(seed)/condition
        inputs[f'{seed}/{condition}']={}
        for filename in ['input.csv','true_counts.csv','shepherd_input.txt']:
            assert sha(ad/filename)==sha(cd/filename), (seed,condition,filename)
            inputs[f'{seed}/{condition}'][filename]=sha(ad/filename)
        truths=set(pd.read_csv(ad/'true_counts.csv').iloc[:,0])
        observed=pd.read_csv(ad/'input.csv')
        absent=len(truths-set(observed.barcode))
        expected_reads=int(observed.counts.sum())
        for idx,row in a[(a.seed==seed)&(a.condition==condition)].iterrows():
            f=ad/f'{row.variant}_{row.method}_{row["repeat"]}.csv'
            result=pd.read_csv(f)
            assert sha(f)==row.output_sha256
            assert row.fn==len(truths-set(result.central_barcode))
            assert row.fp==len(set(result.central_barcode)-truths)
            assert int(result.sum_counts.sum())==expected_reads==row.output_reads
            assert row.absent_truth==absent
        for idx,row in c[(c.seed==seed)&(c.condition==condition)].iterrows():
            f=cd/f'fresh_{row.method}.csv'
            result=pd.read_csv(f)
            assert row.fn==len(truths-set(result.central_barcode))
            assert row.fp==len(set(result.central_barcode)-truths)
            assert int(result.sum_counts.sum())==row.reads
            competitor_outputs[f'{seed}/{condition}/{row.method}']=sha(f)
            c.loc[idx,'input_reads']=expected_reads
            c.loc[idx,'absent_truth']=absent
            c.loc[idx,'fn_present']=row.fn-absent
            c.loc[idx,'output_sha256']=sha(f)
            c.loc[idx,'n_centroids']=len(result)
    c['output_reads']=c.pop('reads')
    c['repeat']=0
    c['wall_s']=c.algorithm_s
    c['f1']=2*(cfg['n_barcodes']-c.fn)/(2*(cfg['n_barcodes']-c.fn)+c.fn+c.fp)
    c['build_id']=''
    combined=pd.concat([a,c],ignore_index=True).sort_values(['seed','condition','variant','method','repeat'])
    assert not combined.duplicated(['seed','condition','variant','method','repeat']).any()
    expected_f1=2*(cfg['n_barcodes']-combined.fn)/(2*(cfg['n_barcodes']-combined.fn)+combined.fn+combined.fp)
    assert ((combined.f1-expected_f1).abs()<1e-12).all()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    combined.to_csv(args.output_dir/'summary_2026-09-08.csv',index=False)
    aggregate=combined.groupby(['condition','variant','method'],as_index=False).agg(
        seeds=('seed','nunique'),mean_fn=('fn','mean'),mean_fp=('fp','mean'),
        mean_ws=('ws','mean'),mean_f1=('f1','mean'),min_f1=('f1','min'),max_f1=('f1','max'),
        median_algorithm_s=('algorithm_s','median'),median_wall_s=('wall_s','median'))
    aggregate.to_csv(args.output_dir/'aggregate_2026-09-08.csv',index=False)
    manifest.update(date='2026-09-08',input_sha256=inputs,competitor_output_sha256=competitor_outputs,
        validation='All 96 rows checked for unique keys, direct set-based FN/FP, input identity, F1 arithmetic, output hashes and barbac read conservation.',
        timing_caveat='Accuracy campaigns ran concurrently; their timings are exploratory. See timing_2026-09-08.csv for separately repeated barbac timings after both accuracy campaigns finished.',
        competitors_measured=['shepherd','starcode'])
    (args.output_dir/'versions_2026-09-08.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(aggregate.to_string(index=False))


if __name__=='__main__': main()
