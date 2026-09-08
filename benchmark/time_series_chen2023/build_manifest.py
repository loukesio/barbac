"""Join ENA runs to the author's primer metadata; select the hBFA1 YPD series."""
import csv
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent


def main():
    source = HERE/'generated/sources'
    with open(source/'All_file_primer_info.csv') as handle:
        primers = {r['Filename']: r for r in csv.DictReader(handle)}
    with open(source/'ena_PRJNA912754.tsv') as handle:
        all_runs = list(csv.DictReader(handle, delimiter='\t'))
    rows = []
    for row in all_runs:
        if not re.fullmatch(r'hBFA1_h1_hBFA1-YPD-R[12]-Time(8|16|24|40)', row['sample_alias']):
            continue
        primer = primers[row['sample_alias']+'_R1.fastq.gz']
        _, _, rep, tp = primer['Library'].split('-')
        assert row['library_layout'] == 'PAIRED' and row['library_strategy'] == 'AMPLICON'
        paths = row['fastq_ftp'].split(';')
        md5s = row['fastq_md5'].split(';'); sizes = row['fastq_bytes'].split(';')
        assert len(paths) == len(md5s) == len(sizes) == 2
        out = dict(run=row['run_accession'], sample=primer['Library'],
                   assay='hBFA1', environment='YPD', replicate=rep, generation=int(tp[4:]),
                   ena_sample_alias=row['sample_alias'], sample_accession=row['sample_accession'],
                   offset1=int(primer['R1_bp_to_BC']), offset2=int(float(primer['R2_bp_to_BC'])),
                   inline1=primer['R1_index'], inline2=primer['R2_index'])
        for mate in (1, 2):
            i = next(i for i, p in enumerate(paths) if p.endswith(f'_{mate}.fastq.gz'))
            out.update({f'r{mate}_url': 'https://'+paths[i], f'r{mate}_md5': md5s[i],
                        f'r{mate}_bytes': int(sizes[i])})
        rows.append(out)
    rows.sort(key=lambda r: (r['replicate'], r['generation']))
    assert len(rows) == 8 and len({r['run'] for r in rows}) == 8
    assert {(r['replicate'], r['generation']) for r in rows} == {
        (rep, t) for rep in ('R1', 'R2') for t in (8, 16, 24, 40)}
    with open(HERE/'samples.tsv', 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    with open(source/'hBFA1_all_freqs_tidy.csv') as handle:
        published = [r for r in csv.DictReader(handle) if r['Test_Environment'] == 'YPD']
    assert len({r['Barcode'] for r in published}) == 2314
    assert len(published) == 18512
    with open(source/'published_YPD.csv', 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(published[0]))
        writer.writeheader(); writer.writerows(published)
    receipt = dict(status='validated', selected_runs=8, paired_fastq_files=16,
        compressed_bytes=sum(r['r1_bytes']+r['r2_bytes'] for r in rows),
        generations=[8, 16, 24, 40], replicates=['R1', 'R2'],
        published_distinct_barcode_pairs=2314, published_rows=len(published),
        note='Published counts are processed reference measurements, not biological ground truth.')
    (HERE/'selection.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
