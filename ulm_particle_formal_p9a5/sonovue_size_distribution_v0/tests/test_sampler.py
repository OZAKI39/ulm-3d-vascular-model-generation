"""Analytic and adversarial cases independent of the real-data smoke thresholds."""
import csv
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sonovue_sampler import HistogramError, SonoVueDistribution, sample_diameters


class SamplerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def histogram(self, rows):
        path = self.root / "case.csv"
        with path.open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["diameter_center_um", "diameter_low_um", "diameter_high_um", "probability"])
            writer.writerows(rows)
        return path

    def test_uniform_interval_analytic_solution(self):
        d = SonoVueDistribution(self.histogram([[1.05, 1.0, 1.1, 1.0]]))
        u = np.array([0.0, 0.125, 0.25, 0.5, 0.75, 1.0])
        np.testing.assert_allclose(d.inverse_cdf(u), 1.0 + 0.1 * u, rtol=0, atol=3e-16)
        np.testing.assert_allclose(d.cdf([0, 1, 1.025, 1.05, 1.1, 2]), [0, 0, .25, .5, 1, 1], atol=2e-15)
        self.assertAlmostEqual(d.target_mean_um(), 1.05)

    def test_zero_bins_and_exact_plateau_convention(self):
        d = SonoVueDistribution(self.histogram([[1.05, 1.0, 1.1, 0], [1.15, 1.1, 1.2, .25],
                                              [1.25, 1.2, 1.3, 0], [1.35, 1.3, 1.4, .75],
                                              [1.45, 1.4, 1.5, 0]]))
        self.assertEqual(len(d.probability), 5)
        self.assertEqual(d.inverse_cdf(0), 1.1)
        self.assertEqual(d.inverse_cdf(.25), 1.3)
        self.assertEqual(d.inverse_cdf(1), 1.4)
        self.assertEqual(d.cdf(1.25), .25)
        self.assertEqual(d.pdf(1.25), 0)
        u = np.r_[0, np.nextafter(.25, 0), .25, np.nextafter(.25, 1), np.nextafter(1., 0)]
        x = d.inverse_cdf(u)
        self.assertTrue(np.all(((x >= 1.1) & (x < 1.2)) | ((x >= 1.3) & (x < 1.4))))
        np.testing.assert_allclose(d.cdf(x), u, rtol=0, atol=2e-15)
        x = d.sample_diameters(20000, 17)
        self.assertTrue(np.all(((x >= 1.1) & (x < 1.2)) | ((x >= 1.3) & (x < 1.4))))

    def test_gap_has_zero_density_without_widening_bins(self):
        d = SonoVueDistribution(self.histogram([[1.05, 1, 1.1, .5], [1.35, 1.3, 1.4, .5]]))
        self.assertEqual(d.cdf(1.2), .5)
        self.assertEqual(d.pdf(1.2), 0)
        self.assertEqual(d.inverse_cdf(.5), 1.3)
        x = d.sample_diameters(20000, 18)
        self.assertFalse(np.any((x >= 1.1) & (x < 1.3)))

    def test_invalid_histograms_rejected(self):
        invalid = {
            "negative_probability": [[1.05, 1, 1.1, -1]],
            "all_zero": [[1.05, 1, 1.1, 0]],
            "nonfinite_probability": [[1.05, 1, 1.1, float("nan")]],
            "infinite_probability": [[1.05, 1, 1.1, float("inf")]],
            "bad_width": [[1.1, 1, 1.2, 1]],
            "zero_width": [[1, 1, 1, 1]],
            "bad_center": [[1.06, 1, 1.1, 1]],
            "unordered": [[1.15, 1.1, 1.2, .5], [1.05, 1, 1.1, .5]],
            "overlap": [[1.05, 1, 1.1, .5], [1.10, 1.05, 1.15, .5]],
            "unresolvable_positive_tail": [[1.05, 1, 1.1, 1], [1.15, 1.1, 1.2, 1e-30]],
        }
        for name, rows in invalid.items():
            with self.subTest(name=name):
                with self.assertRaises(HistogramError):
                    SonoVueDistribution(self.histogram(rows))

    def test_missing_input_and_missing_column(self):
        with self.assertRaises(FileNotFoundError):
            SonoVueDistribution(self.root / "missing.csv")
        p = self.root / "missing_column.csv"
        p.write_text("diameter_center_um,probability\n1.05,1\n")
        with self.assertRaises(HistogramError):
            SonoVueDistribution(p)

    def test_optional_columns_are_not_used(self):
        p = self.histogram([[1.05, 1, 1.1, .25], [1.15, 1.1, 1.2, .75]])
        a = SonoVueDistribution(p).sample_diameters(1000, 19)
        lines = p.read_text().splitlines()
        p.write_text(lines[0] + ",cdf,bin_origin\n" + "\n".join(x + ",NOT_A_NUMBER,unused" for x in lines[1:]) + "\n")
        np.testing.assert_array_equal(a, SonoVueDistribution(p).sample_diameters(1000, 19))

    def test_no_global_rng_mutation_and_explicit_seed(self):
        d = SonoVueDistribution(self.histogram([[1.05, 1, 1.1, 1]]))
        before = np.random.get_state()
        a = d.sample_diameters(1000, 20)
        after = np.random.get_state()
        self.assertEqual(before[0], after[0])
        np.testing.assert_array_equal(before[1], after[1])
        self.assertEqual(before[2:], after[2:])
        np.testing.assert_array_equal(a, d.sample_diameters(1000, 20))
        self.assertFalse(np.array_equal(a, d.sample_diameters(1000, 21)))
        self.assertEqual(d.sample_diameters(0, 20).shape, (0,))
        for n, seed in [(-1, 1), (1.5, 1), (True, 1), (1, None), (1, -1), (1, True)]:
            with self.subTest(n=n, seed=seed), self.assertRaises(ValueError):
                d.sample_diameters(n, seed)
        for u in [-.01, 1.01, float("nan"), float("inf")]:
            with self.subTest(u=u), self.assertRaises(ValueError):
                d.inverse_cdf(u)

    def test_source_input_hash_and_immutable_arrays(self):
        d = SonoVueDistribution()
        self.assertEqual(d.histogram_sha256, "2c9f783c6169421b06c057e16653878e7385139a179682c0648519e72a9f5198")
        with self.assertRaises(ValueError):
            d.probability[0] = 0
        np.testing.assert_array_equal(sample_diameters(50, 22), d.sample_diameters(50, 22))


if __name__ == "__main__":
    unittest.main(verbosity=2)
