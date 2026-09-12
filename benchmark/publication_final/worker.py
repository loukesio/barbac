"""One fresh process producing the same centroid/member exports for each tool."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def run(method, source, dest, tools, library, length):
    config = json.loads((HERE / 'protocol.json').read_text())
    prep = 0.0
    command = None
    if method in ['hamming', 'lv']:
        command = ['Rscript', str(ROOT / 'benchmark/latest_four_conditions/run_barbac.R'),
                   str(library), str(source / 'input.csv'), str(dest), method, '3',
                   'poisson' if method == 'lv' else 'none']
    elif method == 'shepherd':
        command = [sys.executable, str(tools / 'Shepherd/shepherd_t0.py'),
                   '-f', str(dest / 'input.tsv'), '-l', str(length), '-eps', '3', '-bft', '4']
    elif method.startswith('starcode_'):
        command = [str(tools / 'starcode/starcode'), '-d', '3', '-t', '1', '--print-clusters',
                   '-i', str(dest / 'input.tsv'), '-o', str(dest / 'clusters.tsv')]
        command += ['-s'] if method == 'starcode_sphere' else ['-r', '5']
    elif method == 'bartender':
        start = time.perf_counter()
        inputs = pd.read_csv(source / 'input.csv')
        with (dest / 'expanded.csv').open('w') as stream:
            read_id = 0
            for sequence, count in inputs.itertuples(index=False, name=None):
                # Bounded chunks avoid a huge temporary string for an outlier parent.
                for offset in range(0, int(count), 10000):
                    size = min(10000, int(count) - offset)
                    stream.write(''.join(f'{sequence},{i}\n' for i in range(read_id + 1, read_id + size + 1)))
                    read_id += size
        prep = time.perf_counter() - start
        # Exact arguments assembled by bartender_single_com, with explicit step 1.
        # Calling the native program also propagates a failing exit code reliably.
        command = [str(tools / 'bartender-1.1/bartender_single'), str(dest / 'expanded.csv'),
                   str(dest / 'bartender'), '1', '5.0', '5', '1', '1', '3', '1']
    else:
        raise ValueError(method)
    start = time.perf_counter()
    with (dest / 'tool.log').open('w') as stream:
        subprocess.run(command, cwd=dest, stdout=stream, stderr=subprocess.STDOUT,
                       check=True, timeout=config['timeout_seconds'])
    tool_seconds = time.perf_counter() - start
    core = None
    build = None
    if method in ['hamming', 'lv']:
        for line in (dest / 'tool.log').read_text().splitlines():
            if line.startswith('BARBAC_BUILD_ID='):
                build = line.split('=', 1)[1]
            if line.startswith('BARBAC_ALGO_SECONDS='):
                core = float(line.split('=', 1)[1])
        assert build == config['barbac']['build_id'], build
    else:
        counts = pd.read_csv(source / 'input.csv').set_index('barcode').counts
        if method == 'shepherd':
            centroids = pd.read_csv(dest / 'input_pb_freq.csv')
            centroids.columns = ['central_barcode', 'sum_counts']
            raw = pd.read_csv(dest / 'input_seq_clust.csv')
            own = raw[raw.sequence.isin(centroids.central_barcode)]
            assert own.cluster.is_unique and len(own) == len(centroids)
            roots = own.set_index('cluster').sequence
            members = pd.DataFrame({'member': raw.sequence, 'central_barcode': raw.cluster.map(roots)})
            members = members.dropna(subset=['central_barcode'])
        elif method.startswith('starcode_'):
            raw = pd.read_csv(dest / 'clusters.tsv', sep='\t', header=None,
                              names=['central_barcode', 'sum_counts', 'members'])
            centroids = raw[['central_barcode', 'sum_counts']]
            members = raw.assign(member=raw.members.str.split(',')).explode('member')[['member', 'central_barcode']]
        else:
            raw = pd.read_csv(dest / 'bartender_cluster.csv')
            membership = pd.read_csv(dest / 'bartender_barcode.csv')
            centroids = raw.rename(columns={'Center': 'central_barcode', 'time_point_1': 'sum_counts'})[['central_barcode', 'sum_counts']]
            members = pd.DataFrame({'member': membership['Unique.reads'],
                                   'central_barcode': membership['Cluster.ID'].map(raw.set_index('Cluster.ID').Center)})
        members['member_count'] = members.member.map(counts)
        assert not members.isna().any().any()
        centroids.to_csv(dest / 'centroids.csv', index=False)
        members.to_csv(dest / 'members.csv', index=False)
    (dest / 'worker.json').write_text(json.dumps(dict(command=command, tool_seconds=tool_seconds,
         input_preparation_seconds=prep, barbac_core_seconds=core, build_id=build), indent=2) + '\n')


if __name__ == '__main__':
    method, source, dest, tools, library, length = sys.argv[1:]
    run(method, Path(source), Path(dest), Path(tools), Path(library), int(length))
