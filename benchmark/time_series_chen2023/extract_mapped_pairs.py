"""Join barbac_xtr BAM calls to original paired-read quality and UMI metadata.

Iterating original FASTQ order preserves the publication's first-passing UMI
rule. Barcode calls themselves come from the mapped PEAR consensus, so agreement
with the original-mate parser is audited, not assumed.
"""
import argparse
from collections import Counter
import csv
from itertools import zip_longest
import json
from pathlib import Path
import time

from extract_barcodes import fastq, barcode, sha256, write_counts
from run_sample import checked


def component_calls(path):
    calls = {}
    with open(path) as handle:
        for row in csv.DictReader(handle):
            name, bc = row['read_id'], row['barcode']
            if name in calls:
                raise ValueError(f'Duplicate primary read ID in component calls: {name}')
            if len(bc) != int(row['barcode_length']) or not 24 <= len(bc) <= 28:
                raise ValueError('Invalid observed barcode length in BAM extraction')
            calls[name] = bc
    return calls


def join_pairs(r1, r2, components, output, row):
    output = Path(output)
    if (output/'extraction.json').exists():
        raise FileExistsError('Completed paired extraction already exists')
    output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    diverse = component_calls(Path(components)/'diverse_reads.csv')
    environment = component_calls(Path(components)/'environment_reads.csv')
    calls = {key:(d, environment[key]) for key,d in diverse.items() if key in environment}
    component_counts = dict(diverse=len(diverse), environment=len(environment), paired=len(calls))
    del diverse, environment
    stats = Counter(total_pairs=0, short_reads=0, quality_failed=0,
        no_complete_BAM_barcode_pair=0, umi_duplicates=0, umi_conflicts=0,
        retained_molecules=0, non_acgt_molecules=0, raw_pattern_missing_molecules=0,
        BAM_changed_barcode_pair_molecules=0)
    seen, pairs = {}, Counter()
    inline1, inline2 = Counter(), Counter()
    offset1, offset2 = int(row['offset1']), int(row['offset2'])
    for a,b in zip_longest(fastq(r1), fastq(r2)):
        if a is None or b is None or a[0] != b[0]:
            raise ValueError('FASTQ mate identifiers or record counts differ')
        name,s1,q1 = a; _,s2,q2 = b
        pair = calls.pop(name, None)
        stats['total_pairs'] += 1
        inline1[s1[8:14]] += 1; inline2[s2[8:17]] += 1
        if len(s1) < offset1+26 or len(s2) < offset2+26:
            stats['short_reads'] += 1; continue
        quality = q1[offset1:offset1+26]+q2[offset2:offset2+26]
        if sum(ord(c)-33 for c in quality) < 30*52:
            stats['quality_failed'] += 1; continue
        if pair is None:
            stats['no_complete_BAM_barcode_pair'] += 1; continue
        umi = s1[:8]+s2[:8]
        if umi in seen:
            stats['umi_duplicates'] += 1
            stats['umi_conflicts'] += seen[umi] != pair
            continue
        seen[umi] = pair
        pairs[pair] += 1
        stats['retained_molecules'] += 1
        stats['non_acgt_molecules'] += bool(set(''.join(pair))-set('ACGT'))
        original = barcode(s1, offset1), barcode(s2, offset2)
        if None in original: stats['raw_pattern_missing_molecules'] += 1
        elif original != pair: stats['BAM_changed_barcode_pair_molecules'] += 1
    if calls:
        raise ValueError(f'{len(calls)} BAM read IDs are absent from the input FASTQs')
    if 'archive_read_count' in row and stats['total_pairs'] != int(row['archive_read_count']):
        raise ValueError('Full FASTQ record count differs from selected ENA sample')
    if stats['total_pairs'] == 0:
        raise ValueError('No paired FASTQ records processed')
    assert stats['total_pairs'] == sum(stats[k] for k in ['short_reads','quality_failed',
        'no_complete_BAM_barcode_pair','umi_duplicates','retained_molecules'])
    valid = Counter({p:n for p,n in pairs.items() if not set(''.join(p))-set('ACGT')})
    diverse, environment = Counter(), Counter()
    for (d,e),n in valid.items(): diverse[d] += n; environment[e] += n
    for filename,values in [('barcode_pairs.csv',pairs),('pairs_barbac.csv',valid)]:
        with open(output/filename,'w',newline='') as handle:
            writer=csv.writer(handle); writer.writerow(['diverse_barcode','environment_barcode','counts'])
            writer.writerows((d,e,n) for (d,e),n in sorted(values.items(),key=lambda x:(-x[1],x[0])))
    write_counts(output/'diverse_barbac.csv',diverse)
    write_counts(output/'environment_barbac.csv',environment)
    assert sum(diverse.values()) == sum(environment.values()) == sum(valid.values())
    assert sum(valid.values()) == stats['retained_molecules']-stats['non_acgt_molecules']
    receipt=dict(status='complete',profile='chen2023_mapped_PEAR_flanks_24_28_Q30_first_UMI',
        counts_unit='UMI-deduplicated read pairs',sample=row,stats=dict(stats),
        component_read_counts=component_counts,unique_pairs=len(pairs),
        barbac_input_molecules=sum(valid.values()),diverse_unique=len(diverse),environment_unique=len(environment),
        diverse_lengths=dict(Counter({length:sum(n for bc,n in diverse.items() if len(bc)==length) for length in range(24,29)})),
        inline1_top=inline1.most_common(5),inline2_top=inline2.most_common(5),
        seconds=time.perf_counter()-started,quality_cutoff=30,offset1=offset1,offset2=offset2,
        bam_extraction=json.loads((Path(components)/'bam_components.json').read_text()),
        inputs={str(p):sha256(p) for p in (r1,r2)},code_sha256=sha256(__file__),
        outputs={p.name:sha256(p) for p in output.glob('*.csv')})
    (output/'extraction.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--index',type=int,required=True)
    args=parser.parse_args()
    with open(Path(__file__).with_name('samples.tsv')) as handle:
        row=list(csv.DictReader(handle,delimiter='\t'))[args.index]
    raw=args.work_dir/'raw'/row['run']
    files=[raw/f'{row["run"]}_{m}.fastq.gz' for m in (1,2)]
    if not all(checked(p,row,m) for m,p in enumerate(files,1)):
        raise ValueError('Input FASTQ MD5/size verification failed')
    receipt=join_pairs(*files,args.work_dir/'bam_components'/row['sample'],
        args.work_dir/'extracted'/row['sample'],row)
    print(json.dumps(receipt['stats'],indent=2))


if __name__ == '__main__': main()
