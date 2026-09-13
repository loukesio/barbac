"""Run only the separately registered Shepherd sensitivity configuration."""
import argparse
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import pandas as pd

HERE = Path(__file__).resolve().parent
MAIN = Path('/Users/theodosiou/Documents/Projects/test_barbac')
FINAL = MAIN / '.codex/publication-final-2026-09-12/source/benchmark/publication_final'
sys.path.insert(0, str(FINAL))
from common import ENVIRONMENT, sha, save, load, verify_outputs
from run import check_freeze as check_original_freeze
from metrics import score


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def fingerprint():
    cfg = load(HERE / 'protocol.json')
    paths = [HERE / name for name in ['protocol.json', 'run.py', 'worker.py']]
    paths += [Path(cfg['shepherd_script'])]
    paths += [FINAL / name for name in ['freeze.json', 'results.json', 'timing_results.json',
                                      'publication_results.csv', 'schedule.json', 'metrics.py']]
    return {str(p): sha(p) for p in paths}


def freeze():
    check_original_freeze()
    path = HERE / 'freeze.json'
    if path.exists():
        record = load(path)
        assert record['sha256'] == fingerprint(), 'Supplement source changed'
        return record
    record = dict(frozen_at=now(), sha256=fingerprint(), environment=ENVIRONMENT,
                  expected_cells=121, scope='Post hoc sensitivity; not the original confirmatory configuration')
    save(path, record)
    return record


def execute():
    registered = freeze()
    cfg = load(HERE / 'protocol.json')
    schedule = load(FINAL / 'schedule.json')
    rows = []
    for i, block in enumerate(schedule, 1):
        data = block['data']
        source = Path(data['path'])
        verify_outputs(source, data['sha256'])
        key = 'milos' if data['condition'] == 'milos' else f"{data['seed']}/{data['condition']}"
        dest = HERE / 'generated/results' / key
        signature = hashlib.sha256(json.dumps(dict(freeze=registered, data=data), sort_keys=True).encode()).hexdigest()
        if (dest / 'result.json').exists():
            row = load(dest / 'result.json')
            assert row['signature'] == signature
            verify_outputs(dest, row.get('output_sha256', {}))
        else:
            assert shutil.disk_usage(HERE).free > 8 * 2**30
            dest.mkdir(parents=True, exist_ok=False)
            shutil.copyfile(source / 'input.tsv', dest / 'input.tsv')
            cmd = [sys.executable, str(HERE / 'worker.py'), str(source), str(dest),
                   str(data['nominal_length']), data['condition']]
            print('START', i, '/', len(schedule), key, flush=True)
            before = now()
            begin = time.perf_counter()
            row = dict(condition=data['condition'], seed=data.get('seed'), method='shepherd_documented',
                       started_at=before, signature=signature, command=cmd, dataset_path=str(source))
            try:
                with (dest / 'worker.log').open('w') as stream:
                    subprocess.run(cmd, env={**os.environ, **ENVIRONMENT}, stdout=stream,
                                   stderr=subprocess.STDOUT, check=True, timeout=cfg['timeout_seconds'] + 60)
                elapsed = time.perf_counter() - begin
                measured_end = now()
                inputs, truth, labels = (pd.read_csv(source / n) for n in ['input.csv', 'truth.csv', 'labels.csv'])
                metrics = score(pd.read_csv(dest / 'centroids.csv'), pd.read_csv(dest / 'members.csv'),
                                inputs, truth, labels)
                row.update(status='complete', workflow_seconds=elapsed, finished_at=measured_end,
                           **metrics, **load(dest / 'worker.json'),
                           output_sha256={n: sha(dest/n) for n in ['centroids.csv', 'members.csv', 'worker.json']})
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
                row.update(status='failed', finished_at=now(), elapsed_seconds=time.perf_counter()-begin,
                           error=str(exc), output_sha256={})
            save(dest / 'result.json', row)
        rows.append(row)
        save(HERE / 'results.json', rows)
        save(HERE / 'progress.json', dict(completed=len(rows), total=len(schedule), updated_at=now()))
        print('RESULT', i, key, row['status'], 'F1', row.get('f1_percent'),
              'seconds', row.get('workflow_seconds'), flush=True)
    assert len(rows) == registered['expected_cells']
    freeze()
    save(HERE / 'execution_validation.json', dict(completed_at=now(), cells=len(rows),
         successes=sum(r['status']=='complete' for r in rows), failures=sum(r['status']!='complete' for r in rows),
         original_results_and_frozen_sources_preserved=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['freeze', 'run'])
    args = parser.parse_args()
    (HERE / 'generated').mkdir(exist_ok=True)
    with (HERE / 'generated/campaign.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        freeze() if args.action == 'freeze' else execute()
