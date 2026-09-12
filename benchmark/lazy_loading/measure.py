"""Five alternating baseline/candidate fresh startup pairs; ten compatibility runs."""
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
MAIN = Path('/Users/theodosiou/Documents/Projects/test_barbac')
REFERENCE = MAIN / 'benchmark/frozen_reference_v1'
CANDIDATE = MAIN / '.codex/count-learning-speed-2026-09-12/speed-library'
BASELINE = MAIN / 'app/.runtime/library'
sys.path.insert(0, str(REFERENCE))
from run_suite import check_freeze, check_dataset, config, ENVIRONMENT
from simulate import sha


def save(name, obj):
    (HERE / name).write_text(json.dumps(obj, indent=2) + '\n')


def main():
    check_freeze()
    source = HERE.parents[1]
    frozen = dict(source_sha256={name: sha(source / name) for name in
        ['NAMESPACE', 'R/10_barbac_xtr.R', 'src/clustering.cpp', 'R/11_super_cluster2.R',
         'benchmark/lazy_loading/measure.py', 'benchmark/lazy_loading/startup.R']},
        candidate_library_sha256={name: sha(CANDIDATE / 'barbac' / name)
          for name in ['libs/barbac.so', 'R/barbac.rdb', 'R/barbac.rdx', 'NAMESPACE']},
        repetitions=5, order='baseline/candidate on even pairs; reversed on odd pairs',
        workflow_observations=1, timing_claim='Paired startup only. Workflow compatibility timings are descriptive.')
    save('protocol.json', frozen)
    env = {**os.environ, **ENVIRONMENT}
    observations = []
    for pair in range(5):
        choices = [('baseline', BASELINE), ('candidate', CANDIDATE)]
        if pair % 2: choices.reverse()
        for label, library in choices:
            start = time.perf_counter()
            result = subprocess.run(['Rscript', str(HERE / 'startup.R'), str(library)],
                 env=env, capture_output=True, text=True, check=True, timeout=90)
            elapsed = time.perf_counter() - start
            observations.append(dict(pair=pair, variant=label, wall_seconds=elapsed,
                output=result.stdout, stderr=result.stderr))
            print('STARTUP', pair, label, round(elapsed, 3), flush=True)
    medians = {label: statistics.median(r['wall_seconds'] for r in observations if r['variant'] == label)
               for label in ['baseline', 'candidate']}
    save('startup_results.json', dict(observations=observations, medians=medians,
        paired_savings_seconds=[next(r['wall_seconds'] for r in observations if r['pair']==i and r['variant']=='baseline') -
          next(r['wall_seconds'] for r in observations if r['pair']==i and r['variant']=='candidate') for i in range(5)]))
    rows = []
    for condition in config()['conditions'] + ['milos']:
        data = REFERENCE / 'generated/datasets' / condition
        check_dataset(data)
        for method in ['hamming', 'lv']:
            dest = HERE / 'generated' / condition / method
            dest.mkdir(parents=True, exist_ok=False)
            command = ['Rscript', str(MAIN / 'benchmark/latest_four_conditions/run_barbac.R'),
              str(CANDIDATE), str(data / 'input.csv'), str(dest), method, '3', 'poisson' if method == 'lv' else 'none']
            start = time.perf_counter()
            with (dest / 'tool.log').open('w') as f:
                subprocess.run(command, env=env, stdout=f, stderr=subprocess.STDOUT, check=True, timeout=900)
            elapsed = time.perf_counter() - start
            baseline = json.loads((REFERENCE / 'generated/baseline' / condition / method / 'result.json').read_text())
            digests = {name: sha(dest / name) for name in ['centroids.csv', 'members.csv']}
            assert digests == baseline['output_sha256'], (condition, method)
            rows.append(dict(condition=condition, method=method, canonical_outputs_byte_identical=True,
                R_process_seconds=elapsed, output_sha256=digests, command=command))
            save('compatibility.json', rows)
            print('IDENTICAL', condition, method, round(elapsed, 3), flush=True)
    check_freeze()
    assert frozen['source_sha256'] == {name: sha(source / name) for name in frozen['source_sha256']}
    save('validation.json', dict(status='passed', startup_pairs=5, exact_compatibility_cells=len(rows),
        competitor_runs=0, reserved_final_seeds_evaluated=False, main_preserved=True))


if __name__ == '__main__': main()
