#!/usr/bin/env python3
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import web_perf


class WebPerfTest(unittest.TestCase):
    def run_probe(self, static_ms=20):
        calls = []

        def fetch(base, path, *, gzip_ok=False):
            calls.append((base, path, gzip_ok))
            return {
                "status": 200, "bytes": 1_000_000,
                "total_ms": 3000 if "?" in path else static_ms if path in ("/", "/viewer.js") else 20,
                "ttfb_ms": 15, "first_prefix_ms": 100,
                "prefixes": 10, "prefixes_per_s": 3.33, "mb_per_s": 0.33,
            }

        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory:
            json_out = Path(directory) / "perf.json"
            markdown_out = Path(directory) / "perf.md"
            args = ["web_perf.py", "--base-url", "http://viewer.test:5052",
                    "--json-out", str(json_out), "--markdown-out", str(markdown_out)]
            with patch.object(sys, "argv", args), patch.object(web_perf, "fetch", fetch), \
                    contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                status = web_perf.main()
            report = json.loads(json_out.read_text())
            markdown = markdown_out.read_text()
        return status, report, markdown, calls

    def test_fixed_workload_and_artifacts(self):
        status, report, markdown, calls = self.run_probe()
        self.assertEqual(status, 0)
        self.assertTrue(report["ok"])
        self.assertEqual((report["trace_id"], report["epsilon"], report["delta"]), (1, "0.5", "200"))
        self.assertEqual(report["concurrent"]["workers"], 4)
        self.assertEqual(report["concurrent"]["requests"], 8)
        self.assertEqual(len(calls), 44)
        self.assertTrue(all(base == "http://viewer.test:5052" for base, _, _ in calls))
        stream = "/api/trace/1?epsilon=0.5&delta=200"
        for path in ("/", "/viewer.js", "/api/traces", "/api/trace/1/original"):
            self.assertEqual(sum(p == path for _, p, _ in calls), 6)
        self.assertEqual(sum(p == stream and gz for _, p, gz in calls), 6)
        self.assertEqual(sum(p == stream and not gz for _, p, gz in calls), 14)
        self.assertIn("trace 1 (ε=0.5, δ=200)", markdown)
        self.assertIn("✅", markdown)

    def test_fixed_latency_limit_fails(self):
        status, report, markdown, _ = self.run_probe(static_ms=501)
        self.assertEqual(status, 1)
        self.assertFalse(report["ok"])
        failed = {gate["name"] for gate in report["gates"] if not gate["ok"]}
        self.assertEqual(failed, {"static / median", "static /viewer.js median"})
        self.assertIn("❌", markdown)

    def test_removed_overrides_are_rejected_before_requests(self):
        options = {
            "trace-id": "2", "epsilon": "1", "delta": "100", "runs": "1",
            "workers": "1", "rounds": "1", "max-static-ms": "999",
            "max-api-ms": "999", "max-ttfb-ms": "999", "max-first-prefix-ms": "999",
            "max-stream-ms": "999", "min-stream-mbps": "0.01",
            "min-concurrent-rps": "0.01", "max-concurrent-slowest-ms": "999",
        }
        for option, value in options.items():
            with self.subTest(option=option), \
                    patch.object(sys, "argv", ["web_perf.py", f"--{option}", value]), \
                    patch.object(web_perf, "fetch") as fetch, \
                    contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    web_perf.main()
                self.assertEqual(error.exception.code, 2)
                fetch.assert_not_called()


if __name__ == "__main__":
    unittest.main()
