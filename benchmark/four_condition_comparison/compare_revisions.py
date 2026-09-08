"""Paired, multi-seed comparison of installed barbac revisions and competitors.

Install each revision in a separate R library first. Inputs use sequence-based
ordering and identical simulation settings for every method. The runner records
individual repetitions, process wall time, read conservation, F1 and historical
FN/FP/WS metrics. It checkpoints after every result and verifies reused inputs.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import time
from pathlib import Path

import pandas as pd
from run_benchmark import (
    HERE, REPO_ROOT, CONDITIONS, simulate, evaluate, run_starcode, run_shepherd,
    git_revision, SHEPHERD_DIR, STARCODE_BIN,
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def barbac(library, directory, output, method, tie_break):
    start = time.perf_counter()
    proc = subprocess.run([
        'Rscript', str(HERE / 'run_revision_comparison.R'), str(library),
        str(directory / 'input.csv'), str(output), method, tie_break,
    ], capture_output=True, text=True, check=True)
    wall = time.perf_counter() - start
    metrics = dict(line.split('=', 1) for line in proc.stdout.splitlines()
                   if line.startswith(('BARBAC_ALGO_SECONDS=', 'BARBAC_BUILD_ID=')))
    return wall, float(metrics['BARBAC_ALGO_SECONDS']), metrics['BARBAC_BUILD_ID'], proc.stderr


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate-library', type=Path, required=True)
    parser.add_argument('--baseline-library', type=Path)
    parser.add_argument('--seeds', type=int, nargs='+', default=[42, 43, 44])
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--n-barcodes', type=int, default=10000)
    parser.add_argument('--n-reads', type=int, default=1000000)
    parser.add_argument('--competitors', nargs='*', choices=['shepherd', 'starcode'], default=['shepherd', 'starcode'])
    parser.add_argument('--output-dir', type=Path, default=HERE / 'generated' / 'revision-comparison')
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    if args.repeats < 1 or len(set(args.seeds)) != len(args.seeds):
        parser.error('Use positive repetitions and unique seeds')
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    libraries = [('candidate_sequence', args.candidate_library.resolve(), 'sequence'),
                 ('candidate_support', args.candidate_library.resolve(), 'support')]
    if args.baseline_library:
        libraries.insert(0, ('baseline_sequence', args.baseline_library.resolve(), 'sequence'))
    sources = ['src/clustering.cpp', 'R/11_super_cluster2.R',
               'benchmark/indel_experiment/simulate.py',
               'benchmark/indel_experiment/run_experiment.py',
               'benchmark/four_condition_comparison/compare_revisions.py',
               'benchmark/four_condition_comparison/run_revision_comparison.R']
    config = dict(seeds=args.seeds, repeats=args.repeats, n_barcodes=args.n_barcodes,
                  n_reads=args.n_reads, competitors=args.competitors,
                  libraries=[(name, str(lib), tie) for name, lib, tie in libraries],
                  conditions=CONDITIONS, tie_order='sequence', abundance='lognormal', sigma=1.5)
    manifest_path = root / 'manifest.json'
    manifest = dict(configuration=config, git_revision=git_revision(REPO_ROOT),
                    source_sha256={s: digest(REPO_ROOT / s) for s in sources},
                    shepherd_revision=git_revision(SHEPHERD_DIR),
                    starcode_revision=git_revision(STARCODE_BIN.parent),
                    platform=platform.platform(),
                    installed_package_sha256={name: {str(p.relative_to(lib)): digest(p)
                        for p in sorted((lib/'barbac').rglob('*')) if p.is_file()}
                        for name,lib,tie in libraries},
                    python_packages={n: importlib.metadata.version(n) for n in ['numpy','pandas','rapidfuzz','scipy']},
                    timing='Sequential processes; barbac algorithm includes CSV input and clustering; wall includes R startup. Competitor algorithm time is process wall time.',
                    status='running')
    rows = []
    if manifest_path.exists():
        old = json.loads(manifest_path.read_text())
        if not args.resume:
            raise RuntimeError('Output already exists; use --resume or a new output directory')
        if (old['configuration'] != json.loads(json.dumps(config)) or
            old['source_sha256'] != manifest['source_sha256'] or
            old['installed_package_sha256'] != manifest['installed_package_sha256']):
            raise RuntimeError('Configuration or source changed; use a new output directory')
        if (root / 'results.csv').exists():
            rows = pd.read_csv(root / 'results.csv').to_dict('records')
    manifest_path.write_text(json.dumps(manifest, indent=2)+'\n')
    completed = {(r['seed'],r['condition'],r['variant'],r['method'],r['repeat']) for r in rows}
    for seed in args.seeds:
        for condition in CONDITIONS:
            directory = root / str(seed) / condition['slug']
            directory.mkdir(parents=True, exist_ok=True)
            sim_config = dict(n_barcodes=args.n_barcodes, n_reads=args.n_reads,
                              barcode_length=len(condition['template']) if condition['template'] else 20,
                              sub_rate=condition['sub_rate'], ins_rate=condition['ins_rate'],
                              del_rate=condition['del_rate'], sigma=1.5, seed=seed,
                              template=condition['template'], abundance='lognormal', tie_order='sequence')
            marker = directory / 'simulation.json'
            if marker.exists():
                saved = json.loads(marker.read_text())
                if saved['configuration'] != sim_config or any(digest(directory/f) != sha for f,sha in saved['sha256'].items()):
                    raise RuntimeError(f'Input configuration or checksum mismatch: {directory}')
            else:
                simulate(directory, **sim_config)
                marker.write_text(json.dumps(dict(configuration=sim_config,
                    sha256={f: digest(directory/f) for f in ['input.csv','true_counts.csv','shepherd_input.txt']}),indent=2)+'\n')
            input_data = pd.read_csv(directory/'input.csv')
            truth = pd.read_csv(directory/'true_counts.csv')
            input_reads = int(input_data.iloc[:,1].sum())
            absent = len(set(truth.iloc[:,0])-set(input_data.iloc[:,0]))
            assert input_reads == int(truth.iloc[:,1].sum()) == args.n_reads
            tasks=[]
            for repeat in range(args.repeats):
                # Reverse order on alternate repeats to balance cache/order effects.
                for variant, library, tie in libraries[::1 if repeat%2==0 else -1]:
                    for method in ['hamming','lv']:
                        tasks.append((variant, method, repeat, library, tie))
            tasks += [(name,name,0,None,None) for name in args.competitors]
            for variant, method, repeat, library, tie in tasks:
                key=(seed,condition['slug'],variant,method,repeat)
                if key in completed: continue
                output=directory/f'{variant}_{method}_{repeat}.csv'
                build=''; warnings=''
                if library:
                    wall, algorithm, build, warnings = barbac(library,directory,output,method,tie)
                elif method == 'starcode':
                    wall=algorithm=run_starcode(directory/'shepherd_input.txt',output)
                else:
                    centroids,wall=run_shepherd(directory/'shepherd_input.txt',directory,sim_config['barcode_length'])
                    algorithm=wall
                    data=pd.read_csv(centroids)
                    data.columns=['central_barcode','sum_counts']; data.to_csv(output,index=False)
                output_data=pd.read_csv(output)
                assert not output_data.central_barcode.duplicated().any()
                reads=int(output_data.sum_counts.sum())
                if library and reads != input_reads:
                    raise RuntimeError(f'barbac lost reads: {key}')
                output_hash=digest(output)
                previous=next((r for r in rows if r['seed']==seed and r['condition']==condition['slug']
                               and r['variant']==variant and r['method']==method), None)
                if previous:
                    if previous['output_sha256'] != output_hash:
                        raise RuntimeError(f'Non-deterministic repeated output: {key}')
                    fn,fp,ws=previous['fn'],previous['fp'],previous['ws']
                else:
                    scored=evaluate(variant,output,directory/'true_counts.csv',wall,algo_time_s=algorithm)
                    fn,fp,ws=scored.fn,scored.fp,scored.ws
                tp=len(truth)-fn
                row=dict(seed=seed,condition=condition['slug'],variant=variant,method=method,repeat=repeat,
                         fn=fn,fp=fp,ws=ws,f1=2*tp/(2*tp+fn+fp),
                         absent_truth=absent,fn_present=fn-absent,n_centroids=len(output_data),
                         input_reads=input_reads,output_reads=reads,algorithm_s=algorithm,wall_s=wall,
                         build_id=build,output_sha256=output_hash)
                rows.append(row)
                pd.DataFrame(rows).to_csv(root/'results.csv',index=False)
                if warnings: output.with_suffix('.log').write_text(warnings)
                print(json.dumps(row),flush=True)
    manifest['status']='complete'
    manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')


if __name__ == '__main__':
    main()
