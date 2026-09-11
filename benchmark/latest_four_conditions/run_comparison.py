"""Current barbac at distances 2/3 on the four original categories, three seeds.

Revalidates saved Shepherd/Starcode sphere outputs on byte-identical inputs.
No parent labels exist for these inputs: scores describe centroid recovery.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import time

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
BASE = ROOT / 'benchmark/four_condition_comparison'
CONDITIONS = ['random_substitutions', 'random_low_indels',
              'anchored_substitutions', 'anchored_low_indels']
BUILD = 'barbac-2026-09-08-bounded-indels-v13'


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def score(centroids, inputs, truth):
    assert list(centroids.columns) == ['central_barcode', 'sum_counts']
    assert centroids.central_barcode.is_unique and not centroids.isna().any().any()
    assert (centroids.sum_counts > 0).all() and (centroids.sum_counts % 1 == 0).all()
    true_set, found = set(truth.iloc[:, 0]), set(centroids.central_barcode)
    tp = len(true_set & found)
    fn, fp = len(true_set - found), len(found - true_set)
    absent = true_set - set(inputs.barcode)
    recovered = centroids.set_index('central_barcode').sum_counts.reindex(truth.iloc[:, 0], fill_value=0)
    return dict(tp=tp, fn=fn, fp=fp, f1=2*tp/(2*tp+fn+fp),
                precision=tp/len(found), recall=tp/len(true_set),
                absent_truth=len(absent), observed_fn=len((true_set-found)-absent),
                pearson_r=float(pd.Series(recovered.to_numpy()).corr(pd.Series(truth.iloc[:, 1].to_numpy()))),
                centroids=len(found), input_reads=int(inputs.counts.sum()),
                output_reads=int(centroids.sum_counts.sum()),
                unassigned_reads=int(inputs.counts.sum()-centroids.sum_counts.sum()))


def validate_members(directory, inputs, centroids):
    members = pd.read_csv(directory/'members.csv')
    assert members.member.is_unique and not members.isna().any().any()
    assert members.central_barcode.isin(centroids.central_barcode).all()
    pd.testing.assert_series_equal(inputs.set_index('barcode').counts.sort_index(),
                                  members.set_index('member').member_count.sort_index(), check_names=False)
    pd.testing.assert_series_equal(members.groupby('central_barcode').member_count.sum().sort_index(),
                                  centroids.set_index('central_barcode').sum_counts.sort_index(), check_names=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--library', type=Path, default=Path('/private/tmp/barbac-lv-v13'))
    parser.add_argument('--work', type=Path, default=BASE/'generated/latest-four-v13')
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    old = json.loads((ROOT/'benchmark/lv_optimization/manifest_all.json').read_text())
    lib_hashes = {str(p.relative_to(args.library)): sha(p)
                  for p in (args.library/'barbac').rglob('*') if p.is_file()}
    assert lib_hashes == old['libraries'][str(args.library)], 'Use the validated v13 library'
    for path, expected in old['source_sha256'].items():
        assert sha(ROOT/path) == expected, path
    manifest = dict(git_revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                    platform=platform.platform(), library_sha256=lib_hashes,
                    source_sha256={p: sha(ROOT/p) for p in [
                        'src/clustering.cpp', 'R/11_super_cluster2.R',
                        'benchmark/latest_four_conditions/run_barbac.R',
                        'benchmark/latest_four_conditions/run_comparison.py']},
                    seeds=[42,43,44], distances=[2,3], tie_break='support',
                    error_rate=0.005, merge_ratio=20, use_design=False,
                    timing='One serial fresh process per barbac configuration. Core includes CSV reading/order/clustering; workflow adds startup and centroid/member exports. Peer timings are retained earlier serial process observations with centroid output only. No timing confidence intervals.',
                    peers='Shepherd and Starcode sphere, distance 3; saved runs, all input/output hashes revalidated.',
                    parent_labels='Not retained by the original simulator; read assignment accuracy unavailable.',
                    input_sha256={}, output_sha256={}, status='running')
    prior_manifest = args.work/'manifest.json'
    if prior_manifest.exists():
        assert args.resume, 'Use --resume or a fresh output directory'
        saved = json.loads(prior_manifest.read_text())
        assert saved['source_sha256'] == manifest['source_sha256']
        assert saved['library_sha256'] == manifest['library_sha256']
    prior_manifest.write_text(json.dumps(manifest, indent=2)+'\n')
    original = pd.read_csv(BASE/'summary_2026-09-08.csv')
    ablation = pd.read_csv(ROOT/'benchmark/lv_optimization/results_all.csv')
    rows = []
    for seed in [42,43,44]:
        for condition in CONDITIONS:
            source = BASE/'generated/revision-accuracy'/str(seed)/condition
            saved = json.loads((source/'simulation.json').read_text())
            for file, expected in saved['sha256'].items():
                assert sha(source/file) == expected, (source,file)
            manifest['input_sha256'][f'{seed}/{condition}'] = saved
            inputs, truth = pd.read_csv(source/'input.csv'), pd.read_csv(source/'true_counts.csv')
            assert inputs.barcode.is_unique and truth.iloc[:,0].is_unique
            assert len(truth) == 10000 and int(inputs.counts.sum()) == int(truth.iloc[:,1].sum()) == 1000000
            for method in ['shepherd', 'starcode']:
                peer_dir = BASE/'generated/exact-search'/str(seed)/condition
                for file, expected in saved['sha256'].items():
                    assert sha(peer_dir/file) == expected, (peer_dir,file)
                old_row = original[(original.seed==seed) & (original.condition==condition) & (original.variant==method)].iloc[0]
                output = peer_dir/f'fresh_{method}.csv'
                assert sha(output) == old_row.output_sha256, output
                metrics = score(pd.read_csv(output), inputs, truth)
                assert (metrics['fn'], metrics['fp']) == (old_row.fn, old_row.fp)
                assert abs(metrics['f1']-old_row.f1) < 1e-12
                rows.append(dict(seed=seed, condition=condition, method=method, distance=3,
                                 indel_model='not_applicable', origin='retained_peer',
                                 core_seconds=None, workflow_seconds=old_row.wall_s,
                                 output_sha256=sha(output), **metrics))
            # Keep the prior no-Poisson result as a model ablation, with its original time.
            old_row = ablation[(ablation.seed==seed) & (ablation.dataset==condition) & (ablation.variant=='pruned')].iloc[0]
            old_dir = BASE/'generated/lv-v13'/condition/str(seed)/'pruned/support'
            checkpoint = json.loads((old_dir/'result.json').read_text())
            for file, expected in checkpoint['output_sha256'].items():
                assert sha(old_dir/file) == expected
            centroids = pd.read_csv(old_dir/'centroids.csv')
            validate_members(old_dir, inputs, centroids)
            metrics = score(centroids, inputs, truth)
            assert (metrics['fn'], metrics['fp']) == (old_row.fn,old_row.fp)
            rows.append(dict(seed=seed, condition=condition, method='lv', distance=3,
                             indel_model='none', origin='retained_ablation',
                             core_seconds=old_row.core_seconds, workflow_seconds=old_row.process_seconds,
                             output_sha256=sha(old_dir/'centroids.csv'), **metrics))
            tasks = [('hamming',3,'none'), ('lv',3,'poisson'), ('hamming',2,'none'), ('lv',2,'poisson')]
            if seed % 2: tasks.reverse()
            for method, distance, model in tasks:
                key = f'{seed}/{condition}/{method}-d{distance}-{model}'
                dest = args.work/key
                dest.mkdir(parents=True, exist_ok=True)
                checkpoint = dest/'result.json'
                if args.resume and checkpoint.exists():
                    row = json.loads(checkpoint.read_text())
                    for file, expected in row['hashes'].items(): assert sha(dest/file)==expected
                else:
                    print('START',key,flush=True)
                    command = ['Rscript',str(HERE/'run_barbac.R'),str(args.library),
                               str(source/'input.csv'),str(dest),method,str(distance),model]
                    start = time.perf_counter()
                    with (dest/'run.log').open('w') as log:
                        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True,
                                       env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1',
                                            'MKL_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1'})
                    wall = time.perf_counter()-start
                    markers = dict(line.split('=',1) for line in (dest/'run.log').read_text().splitlines()
                                   if line.startswith(('BARBAC_BUILD_ID=','BARBAC_ALGO_SECONDS=')))
                    assert markers['BARBAC_BUILD_ID']==BUILD
                    row = dict(seed=seed, condition=condition, method=method, distance=distance,
                               indel_model=model, origin='fresh_v13', build_id=BUILD,
                               core_seconds=float(markers['BARBAC_ALGO_SECONDS']), workflow_seconds=wall,
                               command=command, hashes={f:sha(dest/f) for f in ['centroids.csv','members.csv']},
                               output_sha256=sha(dest/'centroids.csv'))
                centroids = pd.read_csv(dest/'centroids.csv')
                validate_members(dest, inputs, centroids)
                metrics = score(centroids, inputs, truth)
                if distance==3:
                    if method=='lv':
                        previous = BASE/'generated/lv-v13'/condition/str(seed)/'poisson/support/centroids.csv'
                    else:
                        previous = source/'candidate_support_hamming_0.csv'
                    pd.testing.assert_frame_equal(centroids.sort_values('central_barcode').reset_index(drop=True),
                                                  pd.read_csv(previous).sort_values('central_barcode').reset_index(drop=True))
                    row['prior_centroid_count_equivalence'] = 'passed'
                row.update(metrics)
                checkpoint.write_text(json.dumps(row,indent=2)+'\n')
                manifest['output_sha256'][key]=row['hashes']
                rows.append({k:v for k,v in row.items() if k not in ['hashes','command']})
                pd.DataFrame(rows).to_csv(HERE/'results.csv',index=False)
                print(f"DONE {key}: FN={row['fn']} FP={row['fp']} F1={row['f1']:.6f} core={row['core_seconds']:.3f}s workflow={row['workflow_seconds']:.3f}s",flush=True)
    frame = pd.DataFrame(rows)
    assert len(frame)==84
    frame.to_csv(HERE/'results.csv',index=False)
    aggregate = frame.groupby(['condition','method','distance','indel_model','origin']).agg(
        seeds=('seed','nunique'), fn=('fn','mean'), fp=('fp','mean'), f1=('f1','mean'),
        f1_min=('f1','min'), f1_max=('f1','max'), core_seconds=('core_seconds','mean'),
        workflow_seconds=('workflow_seconds','mean'), unassigned_reads=('unassigned_reads','mean'),
        absent_truth=('absent_truth','mean'), observed_fn=('observed_fn','mean')).reset_index()
    assert (aggregate.seeds==3).all()
    aggregate.to_csv(HERE/'aggregate.csv',index=False)
    manifest['status']='complete'
    manifest['validated_fresh_runs']=48
    manifest['validated_retained_runs']=36
    prior_manifest.write_text(json.dumps(manifest,indent=2)+'\n')
    (HERE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(aggregate.to_string(index=False),flush=True)


if __name__=='__main__': main()
