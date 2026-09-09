"""Check that BAM order cannot change quality filtering or UMI deduplication."""
import csv
import json
from pathlib import Path
import tempfile
import unittest

from extract_mapped_pairs import join_pairs, component_calls
from test_extraction import sequence


class MappedPairsTests(unittest.TestCase):
    def test_original_order_controls_UMIs_and_indels_survive(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder); bc='ACGT'*6+'AC'; indel=bc[:10]+'T'+bc[10:]
            rows=[(bc,'AAAAAAAA','!'),(indel,'AAAAAAAA','I'),
                  (bc,'AAAAAAAA','I'),(bc,'CCCCCCCC','I'),(bc,'GGGGGGGG','I')]
            files=[]
            for mate,offset in [(1,63),(2,49)]:
                f=p/f'r{mate}.fastq'; files.append(f)
                with f.open('w') as h:
                    for i,(b,u,q) in enumerate(rows):
                        s=sequence(b,offset,u)
                        h.write(f'@read{i}/{mate}\n{s}\n+\n{q*len(s)}\n')
            components=p/'components';components.mkdir()
            for filename in ['diverse_reads.csv','environment_reads.csv']:
                with (components/filename).open('w',newline='') as h:
                    writer=csv.writer(h);writer.writerow(['read_id','barcode','barcode_length'])
                    # Deliberately reverse BAM order, and leave read3 unmapped.
                    for i in [4,2,1,0]:writer.writerow([f'read{i}',rows[i][0],len(rows[i][0])])
            (components/'bam_components.json').write_text(json.dumps({'status':'complete'}))
            out=p/'output'
            r=join_pairs(*files,components,out,dict(offset1=63,offset2=49,archive_read_count=5))
            self.assertEqual(r['stats']['quality_failed'],1)
            self.assertEqual(r['stats']['no_complete_BAM_barcode_pair'],1)
            self.assertEqual(r['stats']['umi_duplicates'],1)
            self.assertEqual(r['stats']['umi_conflicts'],1)
            self.assertEqual(r['stats']['retained_molecules'],2)
            with (out/'pairs_barbac.csv').open() as h: calls=list(csv.DictReader(h))
            self.assertEqual({x['diverse_barcode'] for x in calls},{bc,indel})
            # The committed manifest has MD5/size fields, but no archive count.
            without_count=join_pairs(*files,components,p/'without_count',dict(offset1=63,offset2=49))
            self.assertEqual(without_count['stats'],r['stats'])

    def test_duplicate_BAM_identifiers_fail(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'calls.csv';bc='A'*26
            p.write_text(f'read_id,barcode,barcode_length\nr1,{bc},26\nr1,{bc},26\n')
            with self.assertRaisesRegex(ValueError,'Duplicate primary read ID'):
                component_calls(p)


if __name__ == '__main__': unittest.main()
