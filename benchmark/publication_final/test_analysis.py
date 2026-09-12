import unittest
import numpy as np
from analyze import holm, paired_summary


class InferenceTests(unittest.TestCase):
    def test_holm_reference_values_and_order(self):
        np.testing.assert_allclose(holm([.02, .001, .04]), [.04, .003, .04])

    def test_identical_methods_have_no_superiority(self):
        r = paired_summary(np.zeros(30), 1, resamples=1000)
        self.assertFalse(r['superiority_supported'])
        self.assertEqual(r['ties'], 30)
        self.assertEqual(r['one_sided_p'], 1)

    def test_strong_paired_advantage_and_reverse(self):
        differences = np.linspace(.1, .2, 30)
        a = paired_summary(differences, 2, resamples=1000)
        b = paired_summary(-differences, 2, resamples=1000)
        self.assertTrue(a['superiority_supported'])
        self.assertFalse(b['superiority_supported'])
        self.assertEqual(a['wins'], 30)
        self.assertEqual(b['losses'], 30)

    def test_missing_pairs_are_not_silently_dropped(self):
        with self.assertRaises(ValueError):
            paired_summary([1, np.nan, 2], 3, resamples=1000)


if __name__ == '__main__':
    unittest.main()
