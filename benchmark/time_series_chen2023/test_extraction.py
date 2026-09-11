"""Regression checks for paired-read integrity, UMI scope, and indel retention."""
import csv
from pathlib import Path
import tempfile
import unittest

from extract_barcodes import barcode, extract


def sequence(bc, offset, umi='ACGTACGT'):
    return umi+'C'*(offset-14)+'GGTACC'+bc+'ATAACT'+'G'*40


class ExtractionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.base = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)

    def files(self, rows):
        paths = [self.base/f'R{i}.fastq' for i in (1, 2)]
        for mate, offset, path in zip((1, 2), (63, 49), paths):
            with open(path, 'w') as handle:
                for i, (bc, umi, quality) in enumerate(rows):
                    s = sequence(bc, offset, umi)
                    handle.write(f'@read{i}/{mate}\n{s}\n+\n{quality*len(s)}\n')
        return paths

    def test_indel_lengths_are_preserved(self):
        for n in (24, 25, 26, 27, 28):
            bc = ('ACGT'*7)[:n]
            self.assertEqual(barcode(sequence(bc, 63), 63), bc)
        self.assertIsNone(barcode(sequence('C'*22, 63), 63))

    def test_umi_scope_and_quality_count_conservation(self):
        bc = 'ACGT'*6+'AC'
        paths = self.files([(bc, 'AAAAAAAA', 'I'), (bc, 'AAAAAAAA', 'I'),
                            ('T'+bc[1:], 'AAAAAAAA', 'I'),
                            (bc, 'CCCCCCCC', '!'), (bc, 'GGGGGGGG', 'I')])
        r = extract(*paths, self.base/'out')
        self.assertEqual(r['stats']['total_pairs'], 5)
        self.assertEqual(r['stats']['umi_duplicates'], 2)
        self.assertEqual(r['stats']['umi_conflicts'], 1)
        self.assertEqual(r['stats']['quality_failed'], 1)
        self.assertEqual(r['barbac_input_molecules'], 2)
        with open(self.base/'out/diverse_barbac.csv') as handle:
            counts = list(csv.DictReader(handle))
        self.assertEqual(counts[0]['counts'], '2')

    def test_ambiguous_bases_retained_in_audit_but_excluded_from_barbac(self):
        paths = self.files([('ACGT'*6+'AN', 'AAAAAAAA', 'I')])
        r = extract(*paths, self.base/'out')
        self.assertEqual(r['stats']['retained_molecules'], 1)
        self.assertEqual(r['stats']['non_acgt_molecules'], 1)
        self.assertEqual(r['barbac_input_molecules'], 0)

    def test_mismatched_identifiers_fail(self):
        paths = self.files([('ACGT'*6+'AC', 'AAAAAAAA', 'I')])
        paths[1].write_text(paths[1].read_text().replace('@read0/', '@different/'))
        with self.assertRaisesRegex(ValueError, 'identifiers differ'):
            extract(*paths, self.base/'out')

    def test_truncated_fastq_fails(self):
        paths = self.files([('ACGT'*6+'AC', 'AAAAAAAA', 'I')])
        paths[1].write_text(paths[1].read_text()[:-10])
        with self.assertRaisesRegex(ValueError, 'Malformed FASTQ'):
            extract(*paths, self.base/'out')

    def test_extra_mate_fails(self):
        paths = self.files([('ACGT'*6+'AC', 'AAAAAAAA', 'I')])
        paths[1].write_text(paths[1].read_text()*2)
        with self.assertRaisesRegex(ValueError, 'different numbers'):
            extract(*paths, self.base/'out')


if __name__ == '__main__':
    unittest.main()
