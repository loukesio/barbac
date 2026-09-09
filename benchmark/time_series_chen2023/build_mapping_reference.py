"""Reconstruct a masked Chen barcode amplicon from fixed pilot-read segments.

This is an empirical candidate, not an author-deposited reference FASTA.
No published barcode identities are used to construct it.
"""
import argparse
from collections import Counter
import json
from pathlib import Path

from extract_barcodes import fastq, REGEX, sha256

HERE = Path(__file__).resolve().parent
LEFT = 'TTAATATGGACTAAAGGAGGCTTTTGTCGACGGATCCGATATCGGTACC'
MIDDLE = 'ATAACTTCGTATAATGTATGCTATACGAAGTTAT'
REVERSE_LEFT = 'TCGAATTCAAGCTTAGATCTGATATCGGTACC'


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def fixed_segments(path, offset, prefix):
    left, middle = Counter(), Counter()
    total = eligible = 0
    for _, sequence, _ in fastq(path):
        total += 1
        hit = REGEX.match(sequence[offset-10:offset+36])
        if hit is None or len(hit.group(2)) != 26:
            continue
        eligible += 1
        start, end = offset-10+hit.start(2), offset-10+hit.end(2)
        left[sequence[prefix:start]] += 1
        middle[sequence[end:end+34]] += 1
    if not eligible:
        raise ValueError('No nominal-length barcode observations in pilot')
    return dict(total_reads=total, nominal_length_observations=eligible,
                modal_left=left.most_common(1)[0],
                modal_middle=middle.most_common(1)[0])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r1', type=Path, default=HERE/'generated/pilot_input/SRR22757105_1.fastq.gz')
    parser.add_argument('--r2', type=Path, default=HERE/'generated/pilot_input/SRR22757105_2.fastq.gz')
    parser.add_argument('--output', type=Path, default=HERE/'reference')
    args = parser.parse_args()
    a, b = fixed_segments(args.r1, 63, 14), fixed_segments(args.r2, 49, 17)
    if (a['modal_left'][0], a['modal_middle'][0], b['modal_left'][0], b['modal_middle'][0]) != (
            LEFT, MIDDLE, REVERSE_LEFT, reverse_complement(MIDDLE)):
        raise ValueError('Pilot fixed segments differ from the reviewed candidate')
    sequence = LEFT+'N'*26+MIDDLE+'N'*26+reverse_complement(REVERSE_LEFT)
    assert len(sequence) == 167
    args.output.mkdir(parents=True, exist_ok=True)
    fasta = args.output/'chen2023_masked_amplicon.fasta'
    fasta.write_text('>chen2023_barcode_amplicon\n'+sequence+'\n')
    receipt = dict(
        status='empirical candidate; extraction boundaries need indel-aware validation',
        author_deposited_reference=False,
        paper='https://doi.org/10.7554/eLife.92899',
        construct='fixed R1 segment + masked BC2 + shared middle + masked reverse BC1 + reverse-complement fixed R2 segment',
        length=167, reference_name='chen2023_barcode_amplicon',
        coordinates='one-based inclusive',
        barcode_regions=[dict(name='BC2', start=50, end=75, output_orientation='reference'),
                         dict(name='BC1', start=110, end=135, output_orientation='reverse complement of reference')],
        masked_bases=52, barcode_identities_used=False,
        evidence='First 50,000 pairs of SRR22757105; modal fixed segments among nominal 26-base barcode observations',
        r1=a, r2=b,
        inputs={p.name:sha256(p) for p in (args.r1, args.r2)},
        fasta_sha256=sha256(fasta), builder_sha256=sha256(__file__))
    (args.output/'provenance.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
