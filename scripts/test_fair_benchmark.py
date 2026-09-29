#!/usr/bin/env python3
"""Tests the native timer and curve-file contracts, not implementation text."""
import math
from pathlib import Path
import tempfile
import unittest

from fair_benchmark import core_ms, curve_info, summarize, within_bound


class NativeEvidenceTests(unittest.TestCase):
    def test_native_markers_only(self):
        self.assertEqual(core_ms("Loaded 100 points.\nDP_CORE_MS: 0.0142\n", "dp"), 0.0142)
        self.assertEqual(core_ms("SIMPLIFY_CORE_MS: 0.0020\n", "reference"), 0.002)
        for bad in ("elapsed: 3.0", "DP_CORE_MS: nan", "DP_CORE_MS: inf",
                    "DP_CORE_MS: -1", "DP_CORE_MS: 2\nDP_CORE_MS: 3",
                    "SQUISH_CORE_MS: 0.1"):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                core_ms(bad, "dp")

    def test_curve_count_and_finite_coordinates(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "curve.txt"
            path.write_text("2\n0 0\n1.5 -2\n")
            count, identity = curve_info(path)
            self.assertEqual(count, 2)
            self.assertEqual(len(identity), 64)
            for bad in ("3\n0 0\n1 1\n", "0\n", "1\nnan 0\n", "1\n1 2 3\n"):
                path.write_text(bad)
                with self.subTest(bad=bad), self.assertRaises(ValueError):
                    curve_info(path)

    def test_sample_std_and_numerical_slack(self):
        result = summarize([1, 2, 3])
        self.assertEqual(result["mean_ms"], 2)
        self.assertEqual(result["std_ms"], 1)  # Sample, not population SD.
        self.assertTrue(within_bound(300.00000000000017, 300))
        self.assertFalse(within_bound(300.001, 300))
        self.assertFalse(within_bound(math.nan, 300))


if __name__ == "__main__":
    unittest.main()
