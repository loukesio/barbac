"""Paired barcode simulations with explicit censoring outside calibrated support.

Unknown repeat contexts are terminal states for further indel recursion. Their
reads are retained, substituted, labelled and counted. This is a finite-support
simulation assumption, not a claim that biological errors stop at run length 14.
"""
from collections import Counter
from pathlib import Path
import csv
import hashlib
import json
import re

import numpy as np

HERE = Path(__file__).resolve().parent
BASES = np.array(list('ACGT'))


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def load_rates(path=None):
    path = path or HERE / 'homopolymer_rates.csv'
    with Path(path).open() as stream:
        return {int(r['homopolymer_length']):
                (float(r['single_insertion_rate']), float(r['single_deletion_rate']))
                for r in csv.DictReader(stream)}


def embed(identity, template):
    assert template.count('N') == len(identity)
    chars = iter(identity)
    return ''.join(next(chars) if base == 'N' else base for base in template)


def homopolymers(sequence):
    return [(m.start(), m.end(), m.group()[0]) for m in re.finditer(r'(A+|C+|G+|T+)', sequence)
            if m.end() - m.start() >= 5]


def slip(sequence, count, rng, rates):
    """Published-style recursive Poisson allocation with current coordinates.

    Each queue node owns disjoint reads. Subtract allocations before spawning
    descendants, so even a Poisson draw above the remaining count conserves reads.
    """
    out = Counter()
    boundary = Counter()
    stack = [(sequence, int(count), 0)]
    stats = Counter()
    while stack:
        seq, remaining, depth = stack.pop()
        if not remaining:
            continue
        if depth > 100:
            raise RuntimeError('Indel recursion exceeded frozen safety bound')
        stats['max_indel_depth'] = max(stats['max_indel_depth'], depth)
        contexts = homopolymers(seq)
        unsupported = [end - start for start, end, _ in contexts if end - start not in rates]
        if unsupported:
            # Retain the entire disjoint queue node. Do not invent an error rate,
            # discard these reads, change the seed, or select another library.
            out[seq] += remaining
            boundary[seq] += remaining
            stats['boundary_reads'] += remaining
            stats['boundary_states'] += 1
            stats['max_boundary_repeat_length'] = max(stats['max_boundary_repeat_length'], max(unsupported))
            continue
        for start, end, base in contexts:
            if not remaining:
                break
            length = end - start
            insertion, deletion = rates[length]
            plus = min(remaining, int(rng.poisson(insertion * remaining)))
            remaining -= plus
            if plus:
                stack.append((seq[:end] + base + seq[end:], plus, depth + 1))
                stats['insertion_events'] += plus
            minus = min(remaining, int(rng.poisson(deletion * remaining)))
            remaining -= minus
            if minus:
                stack.append((seq[:end - 1] + seq[end:], minus, depth + 1))
                stats['deletion_events'] += minus
        out[seq] += remaining
    assert sum(out.values()) == count
    return +out, stats, boundary


def substitute(sequence, count, probability, rng):
    out = Counter()
    events = 0
    original = np.array(['ACGT'.index(c) for c in sequence], dtype=np.int8)
    for offset in range(0, int(count), 10000):
        size = min(10000, int(count) - offset)
        changed = rng.random((size, len(sequence))) < probability
        events += int(changed.sum())
        dirty = np.flatnonzero(changed.any(axis=1))
        out[sequence] += size - len(dirty)
        if len(dirty):
            values = np.broadcast_to(original, (len(dirty), len(sequence))).copy()
            positions = np.nonzero(changed[dirty])
            values[positions] = (values[positions] + rng.integers(1, 4, len(positions[0]))) % 4
            out.update(''.join(row) for row in BASES[values])
    assert sum(out.values()) == count
    return +out, events


def make_parents(config, seed):
    rng = np.random.default_rng(np.random.SeedSequence([seed, 0]))
    identities = []
    seen = set()
    while len(identities) < config['n_barcodes']:
        identity = ''.join(rng.choice(BASES, 20))
        if identity not in seen:
            identities.append(identity)
            seen.add(identity)
    distribution = config['abundance']
    components = [rng.exponential(mean, distribution[key]) for key, mean in
                  zip(['ordinary', 'medium', 'high'], distribution['means'])]
    abundance = np.concatenate(components)
    assert len(abundance) == len(identities)
    rng.shuffle(abundance)
    counts = rng.poisson(abundance / abundance.sum() * config['expected_reads'])
    return identities, counts


def generate_condition(config, seed, condition, dest, parents=None):
    dest = Path(dest)
    if dest.exists() and any(dest.iterdir()):
        raise FileExistsError(f'Refusing to replace generated data: {dest}')
    dest.mkdir(parents=True, exist_ok=True)
    identities, counts = parents if parents is not None else make_parents(config, seed)
    design, errors = condition.split('_', 1)
    template = config['designs'][design]
    truth = [embed(identity, template) for identity in identities]
    assert len(set(truth)) == config['n_barcodes']
    mixed = errors == 'mixed'
    rates = load_rates()
    labels = Counter()
    boundary_labels = Counter()
    totals = Counter()
    for idx, (parent, count) in enumerate(zip(truth, counts)):
        rng_indel = np.random.default_rng(np.random.SeedSequence([seed, 1, int(design == 'anchored'), idx]))
        rng_sub = np.random.default_rng(np.random.SeedSequence([seed, 2, int(design == 'anchored'), int(mixed), idx]))
        variants, stats, boundary = slip(parent, int(count), rng_indel, rates) if mixed else (Counter({parent: int(count)}), Counter(), Counter())
        totals['max_indel_depth'] = max(totals['max_indel_depth'], stats.pop('max_indel_depth', 0))
        totals['max_boundary_repeat_length'] = max(totals['max_boundary_repeat_length'], stats.pop('max_boundary_repeat_length', 0))
        totals.update(stats)
        totals['reads_with_net_length_change'] += sum(n for seq, n in variants.items() if len(seq) != len(parent))
        for variant, n in sorted(variants.items()):
            assert boundary.get(variant, 0) in (0, n)
            observed, events = substitute(variant, n, config['substitution_probability'], rng_sub)
            totals['substitution_events'] += events
            for sequence, read_count in observed.items():
                labels[(sequence, parent)] += read_count
                if boundary.get(variant, 0):
                    boundary_labels[(sequence, parent)] += read_count
    observed = Counter()
    parent_totals = Counter()
    for (sequence, parent), count in labels.items():
        observed[sequence] += count
        parent_totals[parent] += count
    assert sum(observed.values()) == sum(map(int, counts))
    assert all(parent_totals[parent] == count for parent, count in zip(truth, counts))
    with (dest / 'truth.csv').open('w') as stream:
        writer = csv.writer(stream)
        writer.writerow(['barcode', 'true_count'])
        writer.writerows(zip(truth, map(int, counts)))
    ordered = sorted(observed.items(), key=lambda item: (-item[1], item[0]))
    with (dest / 'input.csv').open('w') as stream:
        writer = csv.writer(stream)
        writer.writerow(['barcode', 'counts'])
        writer.writerows(ordered)
    with (dest / 'input.tsv').open('w') as stream:
        csv.writer(stream, delimiter='\t').writerows(ordered)
    with (dest / 'labels.csv').open('w') as stream:
        writer = csv.writer(stream)
        writer.writerow(['member', 'true_barcode', 'read_count'])
        writer.writerows((seq, parent, n) for (seq, parent), n in sorted(labels.items()))
    with (dest / 'boundary_labels.csv').open('w') as stream:
        writer = csv.writer(stream)
        writer.writerow(['member', 'true_barcode', 'read_count'])
        writer.writerows((seq, parent, n) for (seq, parent), n in sorted(boundary_labels.items()))
    assert sum(boundary_labels.values()) == totals.get('boundary_reads', 0)
    info = dict(condition=condition, seed=seed, nominal_length=len(template),
                variable_bases=template.count('N'), true_barcodes=len(truth),
                zero_read_truth=sum(counts == 0).item(), input_sequences=len(observed),
                input_reads=sum(observed.values()), event_counts=dict(totals),
                truth_parent_totals_reconcile=True,
                boundary_policy='Terminal state for further indels outside measured repeat lengths; retain reads and substitutions, report censoring explicitly',
                boundary_read_fraction=sum(boundary_labels.values()) / max(1, sum(observed.values())),
                sha256={name: sha(dest / name) for name in ['truth.csv', 'input.csv', 'input.tsv', 'labels.csv', 'boundary_labels.csv']})
    (dest / 'dataset.json').write_text(json.dumps(info, indent=2) + '\n')
    return info
