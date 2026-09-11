"""Reproduce the 50k-pair prefix check and compare with the pinned author parser.

Run fetch_sources.py first. This executes only the reviewed parser's function
and class definitions after verifying its pinned SHA-256, not its command-line
program. NumPy/pandas are required only for this independent reference check.
"""
import argparse
import ast
import csv
import gzip
import json
from pathlib import Path
import re
import subprocess
import tempfile
import zlib

from extract_barcodes import extract, sha256

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE/'generated/pilot_reproduction')
    parser.add_argument('--prefix-cache', type=Path, default=HERE/'generated/pilot_prefix_cache')
    args = parser.parse_args()
    import numpy as np
    import pandas as pd
    source = HERE/'generated/sources/PLT_parse_predemultiplexed.py'
    manifest = json.loads((HERE/'sources.json').read_text())
    expected = next(r['sha256'] for r in manifest['files'] if r['filename'] == source.name)
    if sha256(source) != expected:
        raise ValueError('Author parser checksum mismatch')
    inputs = json.loads((HERE/'pilot_inputs.json').read_text())['inputs']
    args.prefix_cache.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='barbac-pilot-') as directory:
        directory = Path(directory)
        paths = []
        for row in inputs:
            partial = args.prefix_cache/f'prefix{row["mate"]}.gz'
            if not partial.exists():
                subprocess.run(['curl', '--fail', '--location', '--silent', '--show-error',
                    '--retry', '3', '--connect-timeout', '30', '--max-time', '120',
                    '--range', row['compressed_range'], row['url'], '--output', str(partial)], check=True)
            if sha256(partial) != row['prefix_sha256']:
                raise ValueError('FASTQ prefix checksum mismatch')
            lines = zlib.decompressobj(31).decompress(partial.read_bytes()).splitlines(keepends=True)
            if len(lines) < 200000:
                raise ValueError('Prefix contains fewer than 50,000 complete reads')
            path = directory/f'sample_R{row["mate"]}.fastq.gz'
            with gzip.open(path, 'wb') as handle:
                handle.write(b''.join(lines[:200000]))
            paths.append(path)
        result = extract(*paths, args.output, max_pairs=50000)
        parsed = ast.parse(source.read_text())
        definitions = [n for n in parsed.body if isinstance(n, (ast.FunctionDef, ast.ClassDef))
            or (isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'MY_REGEX'
                                                for t in n.targets))]
        scope = dict(gzip=gzip, pd=pd, np=np, re=re, csv=csv, QUALITY_CUTOFF=30,
            read_file_base=str(directory)+'/',
            all_primer_info={'sample_R1.fastq.gz': dict(R1_bp_to_BC=63, R2_bp_to_BC=49)})
        exec(compile(ast.Module(body=definitions, type_ignores=[]), str(source), 'exec'), scope)
        original = scope['BcCounter'](pd.DataFrame({'Library': ['pilot']}))
        original.read_files(['sample_R1.fastq.gz'], 'pilot', True)
        reference = {(v[0], v[1]): v[2] for v in original.bc_dict.values()}
        with open(args.output/'barcode_pairs.csv') as handle:
            actual = {(r['diverse_barcode'], r['environment_barcode']): int(r['counts'])
                      for r in csv.DictReader(handle)}
        if actual != reference:
            raise AssertionError('Extracted pairs differ from the author parser')
        assert original.lib_stats['pilot'][:4] == [50000, 421, 614, 11]
        receipt = dict(status='passed', author_parser_exact_pair_counts_match=True,
            author_parser_sha256=sha256(source), extractor_sha256=sha256(HERE/'extract_barcodes.py'),
            scope='First 50000 pairs of SRR22757105 only', stats=result['stats'])
        (args.output/'author_comparison.json').write_text(json.dumps(receipt, indent=2)+'\n')
        print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
