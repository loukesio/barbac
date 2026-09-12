"""Local campaign paths, fingerprints, and shared immutable execution settings."""
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1]
MAIN = Path('/Users/theodosiou/Documents/Projects/test_barbac')
REFERENCE = MAIN / 'benchmark/frozen_reference_v1'
EARLIER = MAIN / '.codex/paired-indel-evidence-2026-09-12/source/benchmark/paired_indels'
LIBRARY = MAIN / '.codex/count-learning-speed-2026-09-12/speed-library'
TOOLS = Path('/Users/theodosiou/Documents/Projects/Barcodes/barbac-benchmark/tools')
ENVIRONMENT = {key: '1' for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS']}
ENVIRONMENT['PYTHONHASHSEED'] = '0'


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2) + '\n')


def load(path):
    return json.loads(Path(path).read_text())


def preserve_reference():
    # Avoid importing another module named simulate or metrics into the campaign.
    frozen = load(REFERENCE / 'freeze.json')
    for path, digest in frozen['source_sha256'].items():
        assert sha(path) == digest, path
    for file in ['src/clustering.cpp', 'R/11_super_cluster2.R']:
        assert sha(SOURCE / file) == sha(MAIN / file), file
    assert sha(MAIN / 'manuscript/submission_bioinformatics/02_new_results/benchmark/benchmark_table.pdf') == frozen['historical_pdf_sha256']
    return frozen


def verify_outputs(directory, hashes):
    for file, digest in hashes.items():
        assert sha(directory / file) == digest, str(directory / file)
