import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd
from common import HERE, REFERENCE, load
import simulate
from metrics import validate_data


class BoundaryModelTests(unittest.TestCase):
    def test_terminal_state_keeps_all_reads(self):
        seq = 'CG' + 'A' * 14 + 'TC'
        out, stats, censored = simulate.slip(seq, 100, np.random.default_rng(11), simulate.load_rates())
        self.assertEqual(out, {seq: 100})
        self.assertEqual(censored, out)
        self.assertEqual(stats['boundary_reads'], 100)

    def test_measured_insertions_reach_explicit_boundary(self):
        seq = 'CG' + 'A' * 13 + 'TC'
        out, stats, censored = simulate.slip(seq, 1000, np.random.default_rng(12), simulate.load_rates())
        self.assertEqual(sum(out.values()), 1000)
        self.assertGreater(stats['boundary_reads'], 0)
        self.assertEqual(sum(censored.values()), stats['boundary_reads'])
        self.assertTrue(all('A' * 14 in s for s in censored))

    def test_zero_count_outside_support_needs_no_rate(self):
        out, stats, censored = simulate.slip('A' * 14, 0, np.random.default_rng(13), simulate.load_rates())
        self.assertEqual(sum(out.values()), 0)
        self.assertFalse(censored)

    def test_supported_slippage_is_identical_to_frozen_generator(self):
        spec = importlib.util.spec_from_file_location('original_simulator', REFERENCE / 'simulate.py')
        original = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(original)
        rates = simulate.load_rates()
        for seed in range(100):
            seq = 'CT' + 'A' * (5 + seed % 5) + 'GC'
            before, stats_before = original.slip(seq, 100, np.random.default_rng(seed), rates)
            after, stats_after, censored = simulate.slip(seq, 100, np.random.default_rng(seed), rates)
            self.assertEqual(before, after)
            self.assertFalse(censored)
            self.assertEqual(dict(stats_before), dict(stats_after))

    def test_boundary_labels_follow_substitution_and_reconcile(self):
        cfg = load(REFERENCE / 'protocol.json')
        cfg['n_barcodes'] = 1
        cfg['designs'] = {'random': 'N' * 20}
        seq = 'CG' + 'A' * 14 + 'TCTG'
        with tempfile.TemporaryDirectory() as tmp:
            meta = simulate.generate_condition(cfg, 17, 'random_mixed', Path(tmp) / 'input', ([seq], np.array([1000])))
            data = Path(tmp) / 'input'
            inp, truth, labels, boundary = (pd.read_csv(data / n) for n in ['input.csv', 'truth.csv', 'labels.csv', 'boundary_labels.csv'])
            self.assertEqual(validate_data(inp, truth, labels)['parent_count_absolute_difference'], 0)
            pd.testing.assert_frame_equal(labels, boundary)
            self.assertEqual(meta['boundary_read_fraction'], 1)
            self.assertGreater(meta['event_counts']['substitution_events'], 0)


if __name__ == '__main__':
    unittest.main()
