#!/usr/bin/env python3
"""Behavioral tests for the NDJSON trace stream (GET /api/trace/<id>).

Covers the viewer data diet transport: identity responses carry the binary's
bytes unchanged, and gzip responses (one compressobj over the whole stream)
decompress to exactly those bytes with the right Content-Encoding/Vary.
"""
import gzip
import shutil
import stat
import sys
import tempfile
import unittest
from pathlib import Path

WEB_DIR = Path(__file__).resolve().parent
REPO_ROOT = WEB_DIR.parent
sys.path.insert(0, str(WEB_DIR))

import server  # noqa: E402

NDJSON_LINES = [
    '{"type":"header","v":2,"stream":[[1.0,2.0]]}\n',
    '{"type":"prefix","data":{"p0_idx":0}}\n',
    '{"type":"done","simplified":[[1.0,2.0]]}\n',
]


def _write_stub_simplify(path: Path) -> None:
    body = "import sys\n" + "".join(
        f"sys.stdout.write({line!r})\n" for line in NDJSON_LINES
    )
    path.write_text("#!/usr/bin/env python3\n" + body, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IEXEC)


class TraceStreamTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="trace-stream-"))
        self.data = self.tmp / "data"
        self.build = self.tmp / "build"
        self.data.mkdir()
        self.build.mkdir()
        self.trace_id = 7
        trace_dir = self.data / str(self.trace_id)
        trace_dir.mkdir()
        (trace_dir / "original.txt").write_text(
            "2\n0.0 0.0\n1.0 1.0\n", encoding="utf-8"
        )
        _write_stub_simplify(self.build / "simplify")

        self._saved = {
            "DATA_DIR": server.DATA_DIR,
            "SIMPLIFY_BIN": server.SIMPLIFY_BIN,
        }
        server.DATA_DIR = self.data
        server.SIMPLIFY_BIN = self.build / "simplify"
        self.client = server.app.test_client()

    def tearDown(self):
        for key, value in self._saved.items():
            setattr(server, key, value)
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_identity_stream_matches_binary_output(self):
        resp = self.client.get(f"/api/trace/{self.trace_id}")
        self.assertEqual(resp.status_code, 200)
        self.assertNotIn("Content-Encoding", resp.headers)
        self.assertEqual(resp.data.decode("utf-8"), "".join(NDJSON_LINES))

    def test_gzip_stream_decompresses_to_binary_output(self):
        resp = self.client.get(
            f"/api/trace/{self.trace_id}", headers={"Accept-Encoding": "gzip"}
        )
        self.assertEqual(resp.status_code, 200)
        self.assertEqual(resp.headers.get("Content-Encoding"), "gzip")
        self.assertIn("Accept-Encoding", resp.headers.get("Vary", ""))
        self.assertEqual(
            gzip.decompress(resp.data).decode("utf-8"), "".join(NDJSON_LINES)
        )


if __name__ == "__main__":
    unittest.main()
