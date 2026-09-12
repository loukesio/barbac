"""Three paired full-worker LV timings on the three proposed main benchmarks."""
import datetime
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time

from measure import HERE, MAIN, REFERENCE, CANDIDATE, BASELINE
from run_suite import check_freeze, check_dataset, ENVIRONMENT, TOOLS
from simulate import sha


def save(name, value):
    (HERE / name).write_text(json.dumps(value, indent=2) + '\n')


def main():
    check_freeze()
    conditions = ['random_mixed', 'anchored_mixed', 'milos']
    protocol = dict(frozen_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        conditions=conditions, method='lv', repetitions=3, order='alternate AB/BA by repetition',
        timing='Identical frozen Python worker wall time, including startup, R, input and output conversion; scoring and hashing excluded',
        source_sha256={str(p): sha(p) for p in [Path(__file__).resolve(), REFERENCE / 'worker.py',
            MAIN / 'benchmark/latest_four_conditions/run_barbac.R']},
        library_sha256={label: {name: sha(library / 'barbac' / name) for name in ['NAMESPACE', 'R/barbac.rdb', 'R/barbac.rdx', 'libs/barbac.so']}
            for label, library in [('baseline', BASELINE), ('candidate', CANDIDATE)]},
        rationale='These three scenarios were selected for proposed publication scope after inspecting development outcomes, before these paired measurements.',
        competitor_runs=0, final_seeds_evaluated=False)
    save('workflow_protocol.json', protocol)
    rows = []
    for condition in conditions:
        source = REFERENCE / 'generated/datasets' / condition
        meta = check_dataset(source)
        expected = json.loads((REFERENCE / 'generated/baseline' / condition / 'lv/result.json').read_text())['output_sha256']
        for pair in range(3):
            variants = [('baseline', BASELINE), ('candidate', CANDIDATE)]
            if pair % 2: variants.reverse()
            for variant, library in variants:
                dest = HERE / 'generated/paired_complete' / condition / f'{pair}_{variant}'
                dest.mkdir(parents=True, exist_ok=False)
                command = [sys.executable, str(REFERENCE / 'worker.py'), 'lv', str(source),
                    str(dest), str(TOOLS), str(library), str(meta['nominal_length'])]
                start = time.perf_counter()
                with (dest / 'worker.log').open('w') as stream:
                    subprocess.run(command, env={**os.environ, **ENVIRONMENT}, stdout=stream, stderr=subprocess.STDOUT, check=True, timeout=960)
                elapsed = time.perf_counter() - start
                assert {name: sha(dest / name) for name in expected} == expected
                row = dict(condition=condition, pair=pair, variant=variant, workflow_seconds=elapsed,
                    worker_command=command, outputs_byte_identical=True, **json.loads((dest / 'worker.json').read_text()))
                rows.append(row)
                save('workflow_observations.json', rows)
                print(condition, pair, variant, round(elapsed, 3), flush=True)
    summary = []
    for condition in conditions:
        variants = {variant: [r['workflow_seconds'] for r in rows if r['condition'] == condition and r['variant'] == variant]
                    for variant in ['baseline', 'candidate']}
        summary.append(dict(condition=condition, medians={v: statistics.median(x) for v, x in variants.items()},
            ranges={v: [min(x), max(x)] for v, x in variants.items()},
            paired_savings_seconds=[a - b for a, b in zip(variants['baseline'], variants['candidate'])]))
    check_freeze()
    assert protocol['source_sha256'] == {p: sha(p) for p in protocol['source_sha256']}
    save('workflow_summary.json', summary)
    save('workflow_validation.json', dict(status='passed', runs=len(rows), paired_repetitions_per_condition=3,
        all_centroids_memberships_counts_identical=True, competitor_runs=0, final_seeds_evaluated=False))


if __name__ == '__main__': main()
