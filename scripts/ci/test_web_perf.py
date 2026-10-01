#!/usr/bin/env python3
import contextlib
import gzip
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
import zlib
from unittest.mock import Mock, patch

import web_perf


class WebPerfTest(unittest.TestCase):
    def run_probe(self, static_ms=20, stream_ok=True, gzip_prefix_ms=100, concurrent_error=False):
        calls = []

        def fetch(base, path, *, gzip_ok=False):
            calls.append((base, path, gzip_ok))
            if concurrent_error and web_perf.threading.current_thread() is not web_perf.threading.main_thread():
                raise ConnectionResetError("connection reset")
            return {
                "status": 200, "bytes": 1_000_000,
                "total_ms": 3000 if "?" in path else static_ms if path in ("/", "/viewer.js") else 20,
                "ttfb_ms": 15, "first_prefix_ms": gzip_prefix_ms if gzip_ok else 100,
                "stream_ok": stream_ok,
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

    def test_failed_streams_fail_all_stream_consumers(self):
        status, report, _, _ = self.run_probe(stream_ok=False)
        self.assertEqual(status, 1)
        failed = {gate["name"] for gate in report["gates"] if not gate["ok"]}
        self.assertTrue({"stream unsuccessful count", "stream(gzip) unsuccessful count",
                         "concurrent unsuccessful workload"}.issubset(failed))
        self.assertEqual(report["concurrent"]["successful_requests"], 0)
        self.assertEqual(report["concurrent"]["requests_per_s"], 0)

    def test_gzip_prefix_latency_fails_and_is_reported(self):
        status, report, _, _ = self.run_probe(gzip_prefix_ms=5000)
        self.assertEqual(status, 1)
        failed = {gate["name"] for gate in report["gates"] if not gate["ok"]}
        self.assertEqual(failed, {"stream(gzip) first-prefix median"})
        self.assertEqual(report["stream_gzip"]["first_prefix_ms"], 5000)

    def test_empty_concurrent_workload_fails_with_artifacts(self):
        status, report, _, _ = self.run_probe(concurrent_error=True)
        self.assertEqual(status, 1)
        cc = report["concurrent"]
        self.assertFalse(cc["all_ok"])
        self.assertEqual(cc["requests"], 0)
        self.assertEqual(len(cc["errors"]), 8)
        self.assertEqual(cc["requests_per_s"], 0)

    def test_partial_concurrent_failure_cannot_pass(self):
        def fetch(base, path):
            if web_perf.threading.current_thread().name == "failed-worker":
                raise ConnectionResetError("connection reset")
            return {"status": 200, "stream_ok": True, "bytes": 1000, "total_ms": 3000}

        thread = web_perf.threading.Thread
        names = iter(["failed-worker", "worker-2", "worker-3", "worker-4"])
        with patch.object(web_perf, "fetch", fetch), \
                patch.object(web_perf.threading, "Thread", side_effect=lambda **kw: thread(name=next(names), **kw)):
            cc = web_perf.concurrent("http://viewer.test", "/api/trace/1", 4, 2)
        self.assertFalse(cc["all_ok"])
        self.assertEqual(cc["requests"], 6)
        self.assertEqual(cc["successful_requests"], 6)
        self.assertEqual(len(cc["errors"]), 2)

    def fetch_chunks(self, chunks, *, compressed=False, content_type="application/x-ndjson"):
        clock = [0.0]
        packets = iter(chunks)

        def read1(size):
            now, data = next(packets, (clock[0], b""))
            clock[0] = now
            return data

        response = Mock(status=200)
        response.getheader.side_effect = {"Content-Type": content_type,
                                          "Content-Encoding": "gzip" if compressed else None}.get
        response.read1.side_effect = read1
        conn = Mock()
        conn.getresponse.return_value = response
        with patch.object(web_perf, "_conn", return_value=(conn, "")), \
                patch.object(web_perf.time, "perf_counter", side_effect=lambda: clock[0]):
            result = web_perf.fetch("http://viewer.test", "/api/trace/1", gzip_ok=compressed)
        conn.close.assert_called_once()
        return result

    def test_fetch_decodes_split_records_for_both_encodings(self):
        body = b'{"type":"header"}\n{"output":[], "type" : "prefix"}\n{"type":"done"}\n'
        for compressed in (False, True):
            wire = gzip.compress(body) if compressed else body
            with self.subTest(compressed=compressed):
                result = self.fetch_chunks([(0.01, wire[:10]), (0.1, wire[10:])], compressed=compressed)
                self.assertTrue(result["stream_ok"])
                self.assertEqual(result["prefixes"], 1)
                self.assertEqual(result["first_prefix_ms"], 100)
                self.assertEqual(result["bytes"], len(wire))
                self.assertEqual(result["ttfb_ms"], 10)

    def test_fetch_rejects_failed_and_incomplete_streams(self):
        start = b'{"type":"header"}\n{"type":"prefix"}\n'
        done = b'{"type":"done"}\n'
        bodies = [start, start + b'{"type":"error"}\n',
                  start + done + b'{"type":"error"}\n', start + done + done,
                  start + b'{"type":"unknown"}\n' + done,
                  start + b'not json\n' + done, start + b'[]\n' + done,
                  b'{"type":"prefix"}\n' + done,
                  b'{"type":"header"}\n' + done]
        for body in bodies:
            for compressed in (False, True):
                with self.subTest(body=body, compressed=compressed):
                    wire = gzip.compress(body) if compressed else body
                    result = self.fetch_chunks([(0.1, wire)], compressed=compressed)
                    self.assertFalse(result["stream_ok"])
        result = self.fetch_chunks([(0.1, gzip.compress(start + done)[:-8])], compressed=True)
        self.assertFalse(result["stream_ok"])
        result = self.fetch_chunks([(0.1, start + done)], content_type="application/json")
        self.assertFalse(result["stream_ok"])

    def test_gzip_wrapper_does_not_count_as_first_prefix(self):
        body = b'{"type":"header"}\n{"type":"prefix"}\n{"type":"done"}\n'
        wire = gzip.compress(body)
        result = self.fetch_chunks([(0.01, wire[:10]), (5.0, wire[10:])], compressed=True)
        self.assertTrue(result["stream_ok"])
        self.assertEqual(result["ttfb_ms"], 10)
        self.assertEqual(result["first_prefix_ms"], 5000)

    def test_gzip_prefix_is_measured_before_stream_completion(self):
        compressor = zlib.compressobj(wbits=31)
        first = compressor.compress(b'{"type":"header"}\n{"type":"prefix"}\n')
        first += compressor.flush(zlib.Z_SYNC_FLUSH)
        last = compressor.compress(b'{"type":"done"}\n') + compressor.flush()
        result = self.fetch_chunks([(0.1, first), (5.0, last)], compressed=True)
        self.assertTrue(result["stream_ok"])
        self.assertEqual(result["first_prefix_ms"], 100)
        self.assertEqual(result["total_ms"], 5000)

    def test_transport_error_closes_connection(self):
        conn = Mock()
        conn.getresponse.side_effect = ConnectionResetError("connection reset")
        with patch.object(web_perf, "_conn", return_value=(conn, "")):
            with self.assertRaises(ConnectionResetError):
                web_perf.fetch("http://viewer.test", "/api/trace/1")
        conn.close.assert_called_once()

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
