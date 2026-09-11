"""Run one checked FASTQ download or extraction task, locally or under SLURM."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import time

from extract_barcodes import extract, sha256


def digest(path):
    h = hashlib.md5()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def checked(path, row, mate):
    return (path.is_file() and path.stat().st_size == int(row[f'r{mate}_bytes'])
            and digest(path) == row[f'r{mate}_md5'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['download', 'map', 'extract', 'extract-original'])
    parser.add_argument('--manifest', type=Path, default=Path(__file__).with_name('samples.tsv'))
    parser.add_argument('--index', type=int, required=True, help='Zero-based manifest row')
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--max-pairs', type=int, default=0)
    parser.add_argument('--fastqc', action='store_true')
    args = parser.parse_args()
    with open(args.manifest) as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if not 0 <= args.index < len(rows):
        raise ValueError('Array index outside manifest')
    row = rows[args.index]
    if args.stage in ('map','extract'):
        if args.max_pairs:
            raise ValueError('BAM workflow uses complete samples; choose index 0 for the pilot')
        from mapped_workflow import run
        receipt=run(args.stage,args.work_dir,args.index,row)
        print(json.dumps({'status':receipt['status'],'sample':row['sample'],'stage':args.stage}))
        return
    raw = args.work_dir/'raw'/row['run']; raw.mkdir(parents=True, exist_ok=True)
    paths = [raw/f'{row["run"]}_{i}.fastq.gz' for i in (1, 2)]
    if args.stage == 'download':
        started = time.perf_counter()
        for mate, path in enumerate(paths, 1):
            if checked(path, row, mate):
                continue
            if path.exists():
                raise ValueError(f'Existing FASTQ fails checksum: {path}; inspect or move it before retrying')
            partial = path.with_suffix(path.suffix+'.part')
            subprocess.run(['curl', '--fail', '--location', '--show-error', '--retry', '5',
                '--connect-timeout', '30', '--max-time', '7200',
                row[f'r{mate}_url'], '--output', str(partial)], check=True)
            if not checked(partial, row, mate):
                raise ValueError(f'Download failed size/MD5 validation: {partial}')
            partial.replace(path)
        (raw/'download.json').write_text(json.dumps(dict(sample=row,
            manifest_sha256=sha256(args.manifest), seconds=time.perf_counter()-started,
            files={p.name: {'md5': digest(p), 'bytes': p.stat().st_size} for p in paths}), indent=2)+'\n')
    else:
        for mate, path in enumerate(paths, 1):
            if not checked(path, row, mate):
                raise ValueError(f'FASTQ missing or checksum mismatch: {path}; run download stage first')
        qc = args.work_dir/'qc'/row['sample']
        qc_receipt = None
        if args.fastqc:
            if not shutil.which('fastqc'):
                raise RuntimeError('FastQC not found; activate the environment configured for the cluster')
            qc.mkdir(parents=True, exist_ok=True)
            version = subprocess.check_output(['fastqc', '--version'], text=True).strip()
            started = time.perf_counter()
            subprocess.run(['fastqc', '--threads', '1', '--outdir', str(qc), *map(str, paths)], check=True)
            qc_receipt = dict(version=version, seconds=time.perf_counter()-started)
        scope = f'pilot_{args.max_pairs}' if args.max_pairs else 'full'
        out = args.work_dir/'extracted'/scope/row['sample']
        receipt = extract(*paths, out, offset1=int(row['offset1']), offset2=int(row['offset2']),
                          max_pairs=args.max_pairs)
        receipt.update(sample=row, manifest_sha256=sha256(args.manifest), fastqc=qc_receipt,
                       python=platform.python_version(), platform=platform.platform(),
                       workflow_sha256=sha256(__file__))
        (out/'extraction.json').write_text(json.dumps(receipt, indent=2)+'\n')
        print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
