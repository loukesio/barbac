"""Fresh-process Shepherd measurement; input staging and scoring are excluded."""
import json
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent


def run(source, dest, length, condition):
    cfg = json.loads((HERE / 'protocol.json').read_text())
    command = [sys.executable, cfg['shepherd_script'], '-f', str(dest / 'input.tsv'),
               '-l', str(length), '-eps', '3', '-bft', '-4']
    if condition != 'milos':
        command += ['-e', str(cfg['simulation_substitution_rate'])]
    start = time.perf_counter()
    with (dest / 'tool.log').open('w') as stream:
        subprocess.run(command, cwd=dest, stdout=stream, stderr=subprocess.STDOUT,
                       check=True, timeout=cfg['timeout_seconds'])
    tool_seconds = time.perf_counter() - start
    counts = pd.read_csv(source / 'input.csv').set_index('barcode').counts
    centroids = pd.read_csv(dest / 'input_pb_freq.csv')
    centroids.columns = ['central_barcode', 'sum_counts']
    raw = pd.read_csv(dest / 'input_seq_clust.csv')
    own = raw[raw.sequence.isin(centroids.central_barcode)]
    assert own.cluster.is_unique and len(own) == len(centroids)
    roots = own.set_index('cluster').sequence
    members = pd.DataFrame({'member': raw.sequence, 'central_barcode': raw.cluster.map(roots)})
    members = members.dropna(subset=['central_barcode'])
    members['member_count'] = members.member.map(counts)
    assert not members.isna().any().any()
    centroids.to_csv(dest / 'centroids.csv', index=False)
    members.to_csv(dest / 'members.csv', index=False)
    (dest / 'worker.json').write_text(json.dumps(dict(command=command, tool_seconds=tool_seconds,
        input_preparation_seconds=0), indent=2) + '\n')


if __name__ == '__main__':
    source, dest, length, condition = sys.argv[1:]
    run(Path(source), Path(dest), int(length), condition)
