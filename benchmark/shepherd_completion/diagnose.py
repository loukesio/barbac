"""Explain estimator failures and development recall without tuning any method."""
import collections
import csv
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
MAIN = Path('/Users/theodosiou/Documents/Projects/test_barbac')
FINAL = MAIN / '.codex/publication-final-2026-09-12/source/benchmark/publication_final'
REFERENCE = MAIN / 'benchmark/frozen_reference_v1/generated'


def read(path):
    with path.open() as stream:
        yield from csv.DictReader(stream)


def estimator_diagnosis():
    script = json.loads((HERE / 'protocol.json').read_text())['shepherd_script']
    spec = importlib.util.spec_from_file_location('shepherd_original', script)
    shepherd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(shepherd)
    failures = json.loads((FINAL / 'failure_audit.json').read_text())
    results = []
    for row in failures:
        length = 20 if row['condition'] == 'random_mixed' else 26
        source = FINAL / 'generated/datasets' / str(row['seed']) / row['condition'] / 'input.tsv'
        with source.open() as stream:
            counts = {s:int(n) for s,n in (line.split() for line in stream) if len(s)==length}
        ordered = sorted(counts, key=counts.get, reverse=True)
        # The installed implementation includes index 500: 501 sequences.
        top = ordered[:501]
        upward = downward = equal = 0
        for sequence in top:
            for i, nucleotide in enumerate(sequence):
                for alternative in 'ACGT':
                    if alternative == nucleotide:
                        continue
                    n = counts.get(sequence[:i] + alternative + sequence[i+1:], 0)
                    if n > counts[sequence]:
                        upward += n
                    elif n < counts[sequence]:
                        downward += n
                    else:
                        equal += n
        total = upward + downward + equal
        denominator = sum(counts[s] for s in top) * length
        reproduced = total / (denominator + total)
        original = shepherd.estimate_rho(ordered, counts, length, 500)
        assert abs(reproduced-original) < 1e-14 and original > .1
        results.append(dict(seed=row['seed'], condition=row['condition'], estimated_rate=original,
            generating_substitution_rate=.004, top_sequences=len(top),
            numerator_upward=upward, numerator_downward=downward, numerator_equal=equal,
            upward_fraction=upward/total, denominator=denominator,
            source=str(source), reproduced_installed_estimator=True))
    (HERE/'estimator_diagnosis.json').write_text(json.dumps(results, indent=2)+'\n')


def recall_diagnosis():
    results = []
    for condition in ['random_substitutions','random_mixed','anchored_substitutions','anchored_mixed','milos']:
        source = REFERENCE / 'datasets' / condition
        truth = {r['barcode']:int(r['true_count']) for r in read(source/'truth.csv')}
        counts = {r['barcode']:int(r['counts']) for r in read(source/'input.csv')}
        outputs = REFERENCE / 'baseline' / condition / 'lv'
        found = {r['central_barcode'] for r in read(outputs/'centroids.csv')}
        mapping = {r['member']:r['central_barcode'] for r in read(outputs/'members.csv')}
        missing = set(truth)-found
        origins = collections.defaultdict(list)
        for row in read(source/'labels.csv'):
            if row['true_barcode'] in missing:
                origins[row['true_barcode']].append((row['member'],int(row['read_count'])))
        categories, examples = collections.Counter(), []
        for barcode in sorted(missing):
            category = ('zero_read_truth' if truth[barcode]==0 else
                        'exact_sequence_unobserved' if barcode not in counts else
                        'observed_identity_merged_or_unassigned')
            categories[category] += 1
            if truth[barcode] > 0:
                examples.append(dict(barcode=barcode,true_count=truth[barcode],category=category,
                    observed_exact_count=counts.get(barcode,0),assigned_to=mapping.get(barcode),
                    origins=[dict(sequence=s,reads=n,root=mapping.get(s)) for s,n in origins[barcode]]))
        results.append(dict(condition=condition,fn=len(missing),categories=dict(categories),
            positive_missed_counts=dict(sorted(collections.Counter(e['true_count'] for e in examples).items())),
            positive_missed_multiple_distinct_descendants=sum(len(e['origins'])>=2 for e in examples),
            positive_examples=examples))
    (HERE/'development_recall_audit.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__ == '__main__':
    estimator_diagnosis()
    recall_diagnosis()
    print('Ten failures reproduced; five development-reference recall decompositions saved.')
