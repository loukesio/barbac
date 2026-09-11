"""Extract paired Chen 2023 barcodes into count tables accepted by barbac.

The paper profile follows the published PLT parser: mean Q30 in the nominal
26-base windows, the published flank regex, then first-observation UMI
deduplication within each library. Reads remain in their sequenced orientation.
No clustering or truth-based filtering is performed here.
"""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
from itertools import zip_longest
import json
from pathlib import Path
import re
import time

PATTERN = (
    r'\D*?(GTACC|GGACC|GGTCC|G.TACC|GG.ACC|GGT.CC|GGTA.C|GGTAC.)'
    r'(\D{24,28})'
    r'(.TAACT|A.AACT|AT.ACT|ATA.CT|ATAA.T|ATAAC|AAACT|ATACT|ATAAT)\D*')
REGEX = re.compile(PATTERN)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def fastq(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='ascii') as handle:
        while True:
            title = handle.readline()
            if not title:
                return
            seq, plus, quality = [handle.readline().rstrip('\r\n') for _ in range(3)]
            if not title.startswith('@') or not plus.startswith('+') or not seq or len(seq) != len(quality):
                raise ValueError(f'Malformed FASTQ record: {path}: {title.strip()}')
            if any(ord(c) < 33 or ord(c) > 126 for c in quality):
                raise ValueError('Invalid Phred+33 quality character')
            name = title.split()[0][1:]
            if name.endswith(('/1', '/2')):
                name = name[:-2]
            yield name, seq.upper(), quality


def barcode(seq, offset):
    hit = REGEX.match(seq[offset-10:offset+36])
    return hit.group(2) if hit else None


def write_counts(path, counts):
    with open(path, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['barcode', 'counts', 'barcode_length'])
        writer.writerows((bc, n, len(bc)) for bc, n in sorted(counts.items(), key=lambda x: (-x[1], x[0])))


def extract(r1, r2, output, offset1=63, offset2=49, quality_cutoff=30, max_pairs=0):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if (output/'extraction.json').exists():
        raise FileExistsError(f'Completed extraction exists: {output}; use a new output directory')
    if min(offset1, offset2) < 10 or max_pairs < 0 or not 0 <= quality_cutoff <= 93:
        raise ValueError('Invalid extraction parameters')
    start = time.perf_counter()
    stats = Counter(total_pairs=0, short_reads=0, quality_failed=0, regex_failed=0,
                    umi_duplicates=0, umi_conflicts=0, retained_molecules=0, non_acgt_molecules=0)
    # UMI scope is library, including both barcode components, matching the paper.
    seen, pairs = {}, Counter()
    inline1, inline2, lengths1, lengths2 = Counter(), Counter(), Counter(), Counter()
    records = zip_longest(fastq(r1), fastq(r2))
    while not max_pairs or stats['total_pairs'] < max_pairs:
        a, b = next(records, (None, None))
        if a is None and b is None:
            break
        if a is None or b is None:
            raise ValueError('Paired FASTQs contain different numbers of records')
        if a[0] != b[0]:
            raise ValueError(f'Paired read identifiers differ: {a[0]} / {b[0]}')
        _, s1, q1 = a; _, s2, q2 = b
        stats['total_pairs'] += 1
        lengths1[len(s1)] += 1; lengths2[len(s2)] += 1
        inline1[s1[8:14]] += 1; inline2[s2[8:17]] += 1
        if len(s1) < offset1+26 or len(s2) < offset2+26:
            stats['short_reads'] += 1
            continue
        quality = q1[offset1:offset1+26] + q2[offset2:offset2+26]
        if sum(ord(c)-33 for c in quality) < quality_cutoff*52:
            stats['quality_failed'] += 1
            continue
        d, e = barcode(s1, offset1), barcode(s2, offset2)
        if d is None or e is None:
            stats['regex_failed'] += 1
            continue
        umi = s1[:8]+s2[:8]
        if umi in seen:
            stats['umi_duplicates'] += 1
            stats['umi_conflicts'] += seen[umi] != (d, e)
            continue
        seen[umi] = (d, e)
        pairs[(d, e)] += 1
        stats['retained_molecules'] += 1
        if set(d+e)-set('ACGT'):
            stats['non_acgt_molecules'] += 1
    if not stats['total_pairs']:
        raise ValueError('No read pairs processed')
    assert stats['total_pairs'] == sum(stats[k] for k in
        ['short_reads', 'quality_failed', 'regex_failed', 'umi_duplicates', 'retained_molecules'])
    diverse, environment, valid_pairs = Counter(), Counter(), Counter()
    for (d, e), count in pairs.items():
        if not set(d+e)-set('ACGT'):
            diverse[d] += count; environment[e] += count; valid_pairs[(d, e)] = count
    for filename, values in [('barcode_pairs.csv', pairs), ('pairs_barbac.csv', valid_pairs)]:
        with open(output/filename, 'w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(['diverse_barcode', 'environment_barcode', 'counts'])
            writer.writerows((d, e, n) for (d, e), n in sorted(values.items(), key=lambda x: (-x[1], x[0])))
    write_counts(output/'diverse_barbac.csv', diverse)
    write_counts(output/'environment_barbac.csv', environment)
    assert sum(diverse.values()) == sum(environment.values()) == sum(valid_pairs.values())
    assert sum(valid_pairs.values()) == stats['retained_molecules']-stats['non_acgt_molecules']
    receipt = dict(status='complete', profile='chen2023_published_24_28',
        counts_unit='UMI-deduplicated read pairs (first passing observation per UMI)',
        max_pairs=max_pairs, full_file_scan=max_pairs == 0, offset1=offset1, offset2=offset2,
        offsets_are_zero_based=True, quality_cutoff=quality_cutoff, stats=dict(stats),
        barbac_input_molecules=sum(valid_pairs.values()), unique_pairs=len(pairs),
        diverse_unique=len(diverse), environment_unique=len(environment),
        read1_lengths=dict(lengths1), read2_lengths=dict(lengths2),
        inline1_top=inline1.most_common(5), inline2_top=inline2.most_common(5),
        extraction_seconds=time.perf_counter()-start,
        inputs={str(Path(p).resolve()): sha256(p) for p in (r1, r2)},
        code_sha256=sha256(__file__))
    length_counts = Counter()
    for bc, count in diverse.items():
        length_counts[len(bc)] += count
    receipt['diverse_lengths'] = dict(sorted(length_counts.items()))
    receipt['outputs'] = {p.name: sha256(p) for p in sorted(output.glob('*.csv'))}
    (output/'extraction.json').write_text(json.dumps(receipt, indent=2)+'\n')
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r1', required=True, type=Path)
    parser.add_argument('--r2', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--offset1', type=int, default=63)
    parser.add_argument('--offset2', type=int, default=49)
    parser.add_argument('--quality-cutoff', type=int, default=30)
    parser.add_argument('--max-pairs', type=int, default=0, help='0 scans both full files')
    args = parser.parse_args()
    print(json.dumps(extract(**vars(args)), indent=2))


if __name__ == '__main__':
    main()
