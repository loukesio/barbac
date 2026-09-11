"""Download pinned public metadata and published counts; no FASTQ downloads."""
import hashlib
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent


def main():
    destination = HERE/'generated/sources'
    destination.mkdir(parents=True, exist_ok=True)
    for row in json.loads((HERE/'sources.json').read_text())['files']:
        target = destination/row['filename']
        if target.exists() and hashlib.sha256(target.read_bytes()).hexdigest() == row['sha256']:
            continue
        partial = target.with_suffix(target.suffix+'.part')
        subprocess.run(['curl', '--fail', '--location', '--silent', '--show-error',
                        '--retry', '4', '--connect-timeout', '30', '--max-time', '600',
                        row['url'], '--output', str(partial)], check=True)
        if hashlib.sha256(partial.read_bytes()).hexdigest() != row['sha256']:
            raise ValueError(f'Source changed: {row["filename"]}; review metadata before proceeding')
        partial.replace(target)
        print(target)


if __name__ == '__main__':
    main()
