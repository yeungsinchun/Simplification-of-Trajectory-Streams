#!/usr/bin/env python3
"""Behavioral tests for GET /api/trace/<id>/compare and related routes."""
import stat
import sys
import tempfile
import unittest
from pathlib import Path

WEB_DIR = Path(__file__).resolve().parent
REPO_ROOT = WEB_DIR.parent
sys.path.insert(0, str(WEB_DIR))

import server  # noqa: E402


def _write_executable(path: Path, body: str) -> None:
    path.write_text("#!/usr/bin/env python3\n" + body, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IEXEC)


class CompareApiTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="compare-api-"))
        self.data = self.tmp / "data"
        self.build = self.tmp / "build"
        self.data.mkdir()
        self.build.mkdir()
        self.trace_id = 7
        self.trace_dir = self.data / str(self.trace_id)
        self.trace_dir.mkdir()
        (self.trace_dir / "original.txt").write_text(
            "4\n0.0 0.0\n1.0 1.0\n2.0 0.0\n3.0 1.0\n",
            encoding="utf-8",
        )

        self._saved = {
            "DATA_DIR": server.DATA_DIR,
            "DP_BIN": server.DP_BIN,
            "SQUISH_BIN": server.SQUISH_BIN,
            "DOTS_BIN": server.DOTS_BIN,
            "SIMPLIFY_BIN": server.SIMPLIFY_BIN,
        }
        self._saved_bins = {
            algo: dict(meta) for algo, meta in server.BASELINE_ALGOS.items()
        }

        server.DATA_DIR = self.data
        server.DP_BIN = self.build / "dp"
        server.SQUISH_BIN = self.build / "squish"
        server.DOTS_BIN = self.build / "dots"
        server.SIMPLIFY_BIN = self.build / "simplify"
        server.BASELINE_ALGOS["dp"]["bin"] = server.DP_BIN
        server.BASELINE_ALGOS["squish"]["bin"] = server.SQUISH_BIN
        server.BASELINE_ALGOS["dots"]["bin"] = server.DOTS_BIN

        self.client = server.app.test_client()

    def tearDown(self):
        for key, value in self._saved.items():
            setattr(server, key, value)
        for algo, meta in self._saved_bins.items():
            server.BASELINE_ALGOS[algo] = meta
        try:
            import shutil
            shutil.rmtree(self.tmp, ignore_errors=True)
        except Exception:
            pass

    def _compare(self, **params):
        return self.client.get(f"/api/trace/{self.trace_id}/compare", query_string=params)

    def test_missing_dp_binary_omits_stale_layer(self):
        stale = "2\n9.0 9.0\n8.0 8.0\n"
        (self.trace_dir / "dp_simplified.txt").write_text(stale, encoding="utf-8")

        resp = self._compare(algorithm="dp", run="1", epsilon="0.5")
        self.assertEqual(resp.status_code, 200)
        payload = resp.get_json()
        self.assertIsNotNone(payload.get("baseline_error"))
        self.assertIn("binary not found", payload["baseline_error"])
        self.assertNotIn("dp", payload["layers"])
        self.assertNotIn("baseline", payload["layers"])
        self.assertIsNone(payload["metrics"]["baseline_points"])
        self.assertIsNone(payload["metrics"]["baseline_core_ms"])
        self.assertIn("original", payload["layers"])

    def test_dp_without_core_ms_token_falls_back_to_elapsed(self):
        _write_executable(
            self.build / "dp",
            r"""
import sys, time
time.sleep(0.05)
src, _eps, out = sys.argv[1], sys.argv[2], sys.argv[3]
text = open(src, encoding="utf-8").read().strip().splitlines()
n = int(text[0])
pts = [line.split()[:2] for line in text[1:1 + n]]
keep = [pts[0], pts[-1]]
with open(out, "w", encoding="utf-8") as fh:
    fh.write(f"{len(keep)}\n")
    for x, y in keep:
        fh.write(f"{x} {y}\n")
print("DP simplification complete.")
print(f"Original: {n} points")
print(f"Simplified: {len(keep)} points")
""",
        )
        resp = self._compare(algorithm="dp", run="1", epsilon="0.9")
        self.assertEqual(resp.status_code, 200)
        payload = resp.get_json()
        self.assertIsNone(payload.get("baseline_error"))
        self.assertIn("dp", payload["layers"])
        self.assertEqual(len(payload["layers"]["dp"]), 2)
        self.assertEqual(payload["metrics"]["baseline_points"], 2)
        self.assertIsInstance(payload["metrics"]["baseline_core_ms"], float)
        self.assertGreaterEqual(payload["metrics"]["baseline_core_ms"], 40.0)

    def test_dp_core_ms_token_is_preferred_when_present(self):
        _write_executable(
            self.build / "dp",
            r"""
import sys
src, _eps, out = sys.argv[1], sys.argv[2], sys.argv[3]
open(out, "w", encoding="utf-8").write("2\n0.0 0.0\n3.0 1.0\n")
print("DP_CORE_MS: 12.5")
""",
        )
        resp = self._compare(algorithm="dp", run="1")
        payload = resp.get_json()
        self.assertIsNone(payload.get("baseline_error"))
        self.assertEqual(payload["metrics"]["baseline_core_ms"], 12.5)

    def test_failed_rerun_does_not_keep_previous_polyline(self):
        (self.trace_dir / "dp_simplified.txt").write_text(
            "3\n0.0 0.0\n1.0 1.0\n3.0 1.0\n",
            encoding="utf-8",
        )
        _write_executable(
            self.build / "dp",
            r"""
import sys
print("boom", file=sys.stderr)
sys.exit(2)
""",
        )
        resp = self._compare(algorithm="dp", run="1", epsilon="0.5")
        payload = resp.get_json()
        self.assertTrue(payload.get("baseline_error"))
        self.assertNotIn("dp", payload["layers"])
        self.assertIsNone(payload["metrics"]["baseline_points"])

    def test_run_dots_query_does_not_select_dots(self):
        resp = self._compare(run_dots="1")
        self.assertEqual(resp.status_code, 200)
        payload = resp.get_json()
        self.assertEqual(payload["algorithm"], "none")
        self.assertIsNone(payload.get("baseline_error"))
        self.assertNotIn("dots", payload["layers"])

    def test_layers_endpoint_is_absent(self):
        resp = self.client.get(f"/api/trace/{self.trace_id}/layers")
        self.assertEqual(resp.status_code, 404)

    def test_frechet_baseline_alias_is_absent(self):
        resp = self.client.get(f"/api/trace/{self.trace_id}/frechet/baseline")
        self.assertEqual(resp.status_code, 400)
        payload = resp.get_json()
        self.assertIn("Unknown curve", payload.get("error", ""))


class RealBaselineBinariesTests(unittest.TestCase):
    """Run the compiled DP/SQUISH CLIs through the compare endpoint."""

    @classmethod
    def setUpClass(cls):
        cls.dp = REPO_ROOT / "build" / "dp"
        cls.squish = REPO_ROOT / "build" / "squish"
        if not cls.dp.exists() or not cls.squish.exists():
            raise unittest.SkipTest("compiled dp/squish binaries are not available")

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="compare-real-"))
        self.data = self.tmp / "data"
        self.data.mkdir()
        self.trace_id = 1
        self.trace_dir = self.data / str(self.trace_id)
        self.trace_dir.mkdir()
        original = REPO_ROOT / "data" / "1" / "original.txt"
        self.trace_dir.joinpath("original.txt").write_text(
            original.read_text(encoding="utf-8"),
            encoding="utf-8",
        )

        self._saved = {
            "DATA_DIR": server.DATA_DIR,
            "DP_BIN": server.DP_BIN,
            "SQUISH_BIN": server.SQUISH_BIN,
        }
        self._saved_bins = {
            algo: dict(meta) for algo, meta in server.BASELINE_ALGOS.items()
        }
        server.DATA_DIR = self.data
        server.DP_BIN = self.dp
        server.SQUISH_BIN = self.squish
        server.BASELINE_ALGOS["dp"]["bin"] = self.dp
        server.BASELINE_ALGOS["squish"]["bin"] = self.squish
        self.client = server.app.test_client()

    def tearDown(self):
        for key, value in self._saved.items():
            setattr(server, key, value)
        for algo, meta in self._saved_bins.items():
            server.BASELINE_ALGOS[algo] = meta
        import shutil
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_real_dp_run_returns_polyline_and_elapsed_ms(self):
        resp = self.client.get(
            f"/api/trace/{self.trace_id}/compare",
            query_string={"algorithm": "dp", "run": "1", "epsilon": "0.9"},
        )
        self.assertEqual(resp.status_code, 200)
        payload = resp.get_json()
        self.assertIsNone(payload.get("baseline_error"), payload.get("baseline_error"))
        pts = payload["layers"]["dp"]
        self.assertGreaterEqual(len(pts), 2)
        self.assertLess(len(pts), payload["metrics"]["original_points"])
        self.assertEqual(payload["metrics"]["baseline_points"], len(pts))
        self.assertIsInstance(payload["metrics"]["baseline_core_ms"], float)
        self.assertGreater(payload["metrics"]["baseline_core_ms"], 0.0)
        self.assertEqual(payload["metrics"]["baseline_label"], "DP")

    def test_real_squish_run_returns_polyline_and_elapsed_ms(self):
        resp = self.client.get(
            f"/api/trace/{self.trace_id}/compare",
            query_string={"algorithm": "squish", "run": "1", "ratio": "0.2"},
        )
        self.assertEqual(resp.status_code, 200)
        payload = resp.get_json()
        self.assertIsNone(payload.get("baseline_error"), payload.get("baseline_error"))
        pts = payload["layers"]["squish"]
        self.assertGreaterEqual(len(pts), 2)
        self.assertEqual(payload["metrics"]["baseline_points"], len(pts))
        self.assertIsInstance(payload["metrics"]["baseline_core_ms"], float)
        self.assertGreater(payload["metrics"]["baseline_core_ms"], 0.0)
        self.assertEqual(payload["metrics"]["baseline_label"], "SQUISH")


if __name__ == "__main__":
    unittest.main()
