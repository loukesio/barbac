"""Freeze, generate, and execute the registered final benchmark exactly once."""
import argparse
import datetime
import fcntl
import hashlib
import importlib.metadata
import os
import platform
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time

import pandas as pd
from common import HERE, SOURCE, MAIN, REFERENCE, LIBRARY, TOOLS, ENVIRONMENT, sha, save, load, preserve_reference, verify_outputs
from simulate import generate_condition, make_parents
from metrics import score, validate_data


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def fingerprints():
    files = [HERE / n for n in ['run.py', 'common.py', 'simulate.py', 'metrics.py', 'worker.py', 'analyze.py',
        'protocol.json', 'final_protocol.json', 'sample_size.json', 'homopolymer_rates.csv',
        'test_simulation.py', 'test_analysis.py', 'simulator_validation.json']]
    files += [SOURCE / n for n in ['src/clustering.cpp', 'R/11_super_cluster2.R', 'R/10_barbac_xtr.R', 'NAMESPACE',
                                   'benchmark/latest_four_conditions/run_barbac.R']]
    files += [LIBRARY / 'barbac' / n for n in ['libs/barbac.so', 'R/barbac.rdb', 'R/barbac.rdx', 'NAMESPACE', 'DESCRIPTION']]
    files += [TOOLS / n for n in ['Shepherd/shepherd_t0.py', 'starcode/starcode', 'bartender-1.1/bartender_single']]
    return {str(p): sha(p) for p in files}


def freeze():
    preserve_reference()
    path = HERE / 'freeze.json'
    if path.exists():
        return check_freeze()
    protocol = load(HERE / 'final_protocol.json')
    assert len(set(protocol['final_seeds'])) == protocol['n_independent_libraries_per_design']
    assert not set(protocol['final_seeds']) & set(load(HERE / 'pilot_protocol.json')['seeds'])
    previous = load(SOURCE / 'benchmark/lazy_loading/protocol.json')
    for name, expected in previous['candidate_library_sha256'].items():
        assert sha(LIBRARY / 'barbac' / name) == expected
    for name in ['NAMESPACE', 'R/10_barbac_xtr.R', 'src/clustering.cpp', 'R/11_super_cluster2.R']:
        assert sha(SOURCE / name) == previous['source_sha256'][name]
    record = dict(frozen_at=now(), source_git=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=SOURCE, text=True).strip(),
        sha256=fingerprints(), machine=platform.machine(), platform=platform.platform(), python=sys.version,
        python_packages={name: importlib.metadata.version(name) for name in ['numpy', 'pandas', 'scipy']},
        environment=ENVIRONMENT, expected_generated_inputs=2*len(protocol['final_seeds']),
        expected_timed_cells=(2*len(protocol['final_seeds'])+1)*len(protocol['methods']),
        competitor_runs_are_final_campaign=True,
        reference_receipts={str(p): sha(p) for p in (REFERENCE / 'generated/baseline').glob('*/*/result.json')})
    save(path, record)
    print('FROZEN', record['frozen_at'], record['source_git'], flush=True)
    return record


def check_freeze():
    preserve_reference()
    record = load(HERE / 'freeze.json')
    assert record['sha256'] == fingerprints(), 'Final source changed; retain existing campaign and investigate'
    for path, digest in record['reference_receipts'].items():
        assert sha(path) == digest
    return record


def generate():
    frozen = check_freeze()
    cfg, protocol = load(HERE / 'protocol.json'), load(HERE / 'final_protocol.json')
    records = []
    for seed in protocol['final_seeds']:
        parents = make_parents(cfg, seed)
        for condition in protocol['designs']:
            dest = HERE / 'generated/datasets' / str(seed) / condition
            if (dest / 'dataset.json').exists():
                metadata = load(dest / 'dataset.json')
                verify_outputs(dest, metadata['sha256'])
            else:
                metadata = generate_condition(cfg, seed, condition, dest, parents)
            inputs, truth, labels = (pd.read_csv(dest / n) for n in ['input.csv', 'truth.csv', 'labels.csv'])
            quality = validate_data(inputs, truth, labels)
            assert quality['parent_count_absolute_difference'] == 0
            records.append(dict(**metadata, **quality, path=str(dest)))
            save(HERE / 'datasets.json', records)
            save(HERE / 'progress.json', dict(phase='generation', completed=len(records), total=frozen['expected_generated_inputs'],
                current_seed=seed, current_condition=condition, updated_at=now()))
            print('DATASET', len(records), '/', frozen['expected_generated_inputs'], seed, condition,
                  'boundary reads', metadata['event_counts'].get('boundary_reads', 0), flush=True)
            del inputs, truth, labels
    assert len(records) == frozen['expected_generated_inputs']
    save(HERE / 'generation_validation.json', dict(inputs=len(records), all_reads_and_parent_totals_reconcile=True,
        boundary_reads=sum(r['event_counts'].get('boundary_reads', 0) for r in records),
        inputs_reaching_boundary=sum(r['event_counts'].get('boundary_reads', 0)>0 for r in records),
        total_reads=sum(r['input_reads'] for r in records), all_registered_seeds_retained=True))
    check_freeze()


def execute():
    frozen = check_freeze()
    protocol = load(HERE / 'final_protocol.json')
    datasets = load(HERE / 'datasets.json')
    assert len(datasets) == frozen['expected_generated_inputs']
    rng = random.Random(protocol['execution']['order_seed'])
    rng.shuffle(datasets)
    # A single fixed published dataset is measured in this final session as well.
    milo = REFERENCE / 'generated/datasets/milos'
    datasets.insert(len(datasets)//2, dict(**load(milo/'dataset.json'), path=str(milo)))
    schedule = []
    for data in datasets:
        methods = list(protocol['methods'])
        rng.shuffle(methods)
        schedule.append(dict(data=data, methods=methods))
    save(HERE / 'schedule.json', schedule)
    rows = []
    started = now()
    for block in schedule:
        data = block['data']
        source = Path(data['path'])
        verify_outputs(source, data['sha256'])
        inputs, truth, labels = (pd.read_csv(source / n) for n in ['input.csv', 'truth.csv', 'labels.csv'])
        scope = 'published_reference' if data['condition'] == 'milos' else 'final'
        for method in block['methods']:
            if shutil.disk_usage(HERE).free < protocol['execution']['minimum_free_disk_gib'] * 2**30:
                raise RuntimeError('Disk reserve reached; completed measurements retained')
            key = 'milos' if scope == 'published_reference' else f"{data['seed']}/{data['condition']}"
            dest = HERE / 'generated/results' / key / method
            signature = hashlib.sha256(json.dumps(dict(freeze=frozen, dataset=data, method=method), sort_keys=True).encode()).hexdigest()
            receipt = dest / 'result.json'
            if receipt.exists():
                row = load(receipt)
                assert row['signature'] == signature
                verify_outputs(dest, row.get('output_sha256', {}))
                print('REUSE', key, method, flush=True)
            else:
                dest.mkdir(parents=True, exist_ok=False)
                shutil.copyfile(source / 'input.tsv', dest / 'input.tsv')
                cmd = [sys.executable, str(HERE / 'worker.py'), method, str(source), str(dest), str(TOOLS), str(LIBRARY), str(data['nominal_length'])]
                before = now()
                print('START', len(rows)+1, '/', frozen['expected_timed_cells'], key, method, flush=True)
                begin = time.perf_counter()
                try:
                    with (dest / 'worker.log').open('w') as stream:
                        subprocess.run(cmd, env={**os.environ, **ENVIRONMENT}, stdout=stream, stderr=subprocess.STDOUT, check=True,
                                       timeout=protocol['execution']['timeout_seconds']+60)
                    elapsed = time.perf_counter()-begin
                except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
                    row = dict(status='failed', condition=data['condition'], seed=data.get('seed'), method=method, scope=scope,
                        signature=signature, error=str(exc), command=cmd, started_at=before, finished_at=now(), output_sha256={})
                    save(receipt, row)
                    rows.append(row)
                    save(HERE / 'results.json', rows)
                    continue
                metrics = score(pd.read_csv(dest/'centroids.csv'), pd.read_csv(dest/'members.csv'), inputs, truth, labels)
                output_sha = {n: sha(dest/n) for n in ['centroids.csv','members.csv','worker.json']}
                row = dict(status='complete', condition=data['condition'], seed=data.get('seed'), method=method, scope=scope,
                    signature=signature, started_at=before, finished_at=now(), workflow_seconds=elapsed,
                    worker_command=cmd, **metrics, **load(dest/'worker.json'), output_sha256=output_sha,
                    boundary_read_fraction=data.get('boundary_read_fraction', 0))
                if scope == 'published_reference':
                    old = load(REFERENCE/'generated/baseline/milos'/method/'result.json')
                    row['milo_cached_accuracy_equal'] = all(abs(metrics[k]-old[k]) < 1e-10 for k in metrics)
                    row['milo_cached_mapping_byte_equal'] = all(output_sha[n] == old['output_sha256'][n] for n in ['centroids.csv','members.csv'])
                save(receipt, row)
                if (dest/'expanded.csv').exists():
                    (dest/'expanded.csv').unlink()
            rows.append(row)
            save(HERE/'results.json', rows)
            save(HERE/'progress.json', dict(phase='execution', completed=len(rows), total=frozen['expected_timed_cells'],
                latest_condition=data['condition'], latest_method=method, started_at=started, updated_at=now()))
            print('RESULT', len(rows), key, method, row['status'],
                  'F1', row.get('f1_percent'), 'seconds', row.get('workflow_seconds'), flush=True)
        del inputs, truth, labels
    assert len(rows) == frozen['expected_timed_cells']
    check_freeze()
    save(HERE/'execution_validation.json', dict(cells=len(rows), successful_cells=sum(r['status']=='complete' for r in rows),
        failed_cells=sum(r['status']!='complete' for r in rows), frozen_sources_and_old_caches_preserved=True,
        all_successful_counts_reconcile=True, completed_at=now()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['freeze','generate','run','all'])
    args = parser.parse_args()
    (HERE/'generated').mkdir(exist_ok=True)
    with (HERE/'generated/campaign.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        lock.write(str(os.getpid())); lock.flush()
        if args.action in ['freeze','all']: freeze()
        if args.action in ['generate','all']: generate()
        if args.action in ['run','all']: execute()
