#!/usr/bin/env python3
"""
Latency / throughput probe for the running Flask web viewer (stdlib only).

Complements scripts/ci/validate_web_trace.py (structure) by measuring what a
user of the viewer feels:

  static    GET /, /viewer.js           time to last byte
  api       GET /api/traces, /api/trace/<id>/original
  stream    GET /api/trace/<id>?epsilon&delta   (NDJSON green-path stream)
              ttfb_ms         request sent -> first body byte (header line)
              first_prefix_ms request sent -> first {"type":"prefix"} line
              total_ms        request sent -> last byte
              MB/s, prefixes/s over the whole stream
  concurrent  N parallel stream requests: aggregate requests/s and MB/s,
              and the slowest request's total_ms

Each metric is sampled `--runs` times (first run is a discarded warm-up) and
the median and worst are compared to generous absolute ceilings, which are
tuned to catch order-of-magnitude regressions (buffering the whole stream
before the first byte, serialising concurrent requests, a slow static path)
on noisy shared runners rather than small drifts.

Usage:
  python3 scripts/ci/web_perf.py --base-url http://127.0.0.1:5051
  python3 scripts/ci/web_perf.py --base-url ... --json-out perf.json --markdown-out perf.md

Exit 0 when every gate passes, 1 otherwise (all gates are always evaluated).
"""
from __future__ import annotations

import argparse
import http.client
import json
import statistics
import sys
import threading
import time
import urllib.parse


def _conn(base: str, timeout: float) -> tuple[http.client.HTTPConnection, str]:
    u = urllib.parse.urlsplit(base)
    return http.client.HTTPConnection(u.hostname, u.port or 80, timeout=timeout), u.path.rstrip("/")


def fetch(base: str, path: str, *, timeout: float = 60.0, gzip_ok: bool = False) -> dict:
    """GET path; return timings (ms), size and, for NDJSON, first-prefix time."""
    conn, prefix = _conn(base, timeout)
    headers = {"Accept-Encoding": "gzip" if gzip_ok else "identity"}
    t0 = time.perf_counter()
    conn.request("GET", prefix + path, headers=headers)
    resp = conn.getresponse()
    ttfb = first_prefix = None
    nbytes = prefixes = 0
    tail = b""
    ndjson = "ndjson" in (resp.getheader("Content-Type") or "")
    gz = (resp.getheader("Content-Encoding") or "") == "gzip"
    while True:
        chunk = resp.read1(65536)
        if not chunk:
            break
        now = time.perf_counter()
        if ttfb is None:
            ttfb = now - t0
        nbytes += len(chunk)
        # Timing of the first prefix needs plaintext; only scan identity bodies.
        if ndjson and not gz:
            buf = tail + chunk
            *lines, tail = buf.split(b"\n")
            for ln in lines:
                if ln.startswith(b'{"type":"prefix"') or ln.startswith(b'{"type": "prefix"'):
                    prefixes += 1
                    if first_prefix is None:
                        first_prefix = now - t0
    total = time.perf_counter() - t0
    conn.close()
    out = {"status": resp.status, "bytes": nbytes, "total_ms": total * 1e3,
           "ttfb_ms": (ttfb if ttfb is not None else total) * 1e3}
    if ndjson and not gz:
        out["prefixes"] = prefixes
        if first_prefix is not None:
            out["first_prefix_ms"] = first_prefix * 1e3
        out["prefixes_per_s"] = prefixes / total if total > 0 else 0.0
    out["mb_per_s"] = nbytes / 1e6 / total if total > 0 else 0.0
    return out


def sample(fn, runs: int) -> list[dict]:
    """Run fn() runs+1 times, dropping the first (warm-up: page cache, imports)."""
    results = [fn() for _ in range(runs + 1)]
    return results[1:]


def concurrent(base: str, path: str, workers: int, rounds: int) -> dict:
    """`workers` simultaneous stream requests, repeated `rounds` times."""
    t0 = time.perf_counter()
    results: list[dict] = []
    lock = threading.Lock()

    def worker():
        for _ in range(rounds):
            r = fetch(base, path)
            with lock:
                results.append(r)

    threads = [threading.Thread(target=worker) for _ in range(workers)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    wall = time.perf_counter() - t0
    return {
        "requests": len(results),
        "all_ok": all(r["status"] == 200 for r in results),
        "wall_s": wall,
        "requests_per_s": len(results) / wall,
        "mb_per_s": sum(r["bytes"] for r in results) / 1e6 / wall,
        "max_total_ms": max(r["total_ms"] for r in results),
        "median_total_ms": statistics.median(r["total_ms"] for r in results),
    }


class Gates:
    def __init__(self):
        self.rows: list[tuple[str, float, str, float, bool]] = []

    def le(self, name: str, value: float, limit: float, unit: str = "ms"):
        self.rows.append((name, value, f"≤ {limit:g} {unit}", limit, value <= limit))

    def ge(self, name: str, value: float, limit: float, unit: str = ""):
        self.rows.append((name, value, f"≥ {limit:g} {unit}".strip(), limit, value >= limit))

    @property
    def ok(self) -> bool:
        return all(r[4] for r in self.rows)


def med(rs: list[dict], key: str) -> float:
    return statistics.median(r[key] for r in rs)


def worst(rs: list[dict], key: str) -> float:
    return max(r[key] for r in rs)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-url", default="http://127.0.0.1:5051")
    ap.add_argument("--trace-id", type=int, default=1)
    ap.add_argument("--epsilon", default="0.5")
    ap.add_argument("--delta", default="200")
    ap.add_argument("--runs", type=int, default=5, help="measured runs per metric (plus one warm-up)")
    ap.add_argument("--workers", type=int, default=4, help="parallel stream clients")
    ap.add_argument("--rounds", type=int, default=2, help="requests per parallel client")
    ap.add_argument("--json-out")
    ap.add_argument("--markdown-out")
    # Ceilings are generous (well above a healthy run); see .github/workflows/README.md.
    ap.add_argument("--max-static-ms", type=float, default=500)
    ap.add_argument("--max-api-ms", type=float, default=500)
    ap.add_argument("--max-ttfb-ms", type=float, default=1000)
    ap.add_argument("--max-first-prefix-ms", type=float, default=3000)
    ap.add_argument("--max-stream-ms", type=float, default=10000)
    ap.add_argument("--min-stream-mbps", type=float, default=0.2)
    ap.add_argument("--min-concurrent-rps", type=float, default=0.1)
    ap.add_argument("--max-concurrent-slowest-ms", type=float, default=60000)
    a = ap.parse_args()

    base = a.base_url
    qs = f"epsilon={a.epsilon}&delta={a.delta}"
    stream_path = f"/api/trace/{a.trace_id}?{qs}"
    report: dict = {"trace_id": a.trace_id, "epsilon": a.epsilon, "delta": a.delta}
    g = Gates()

    for name, path in (("static /", "/"), ("static /viewer.js", "/viewer.js"),
                       ("api /api/traces", "/api/traces"),
                       (f"api /api/trace/{a.trace_id}/original", f"/api/trace/{a.trace_id}/original")):
        rs = sample(lambda p=path: fetch(base, p), a.runs)
        limit = a.max_static_ms if name.startswith("static") else a.max_api_ms
        bad = [r["status"] for r in rs if r["status"] != 200]
        g.le(f"{name} status!=200 count", len(bad), 0, "")
        g.le(f"{name} median", med(rs, "total_ms"), limit)
        g.le(f"{name} worst", worst(rs, "total_ms"), limit * 3)
        report[name] = {"median_ms": med(rs, "total_ms"), "worst_ms": worst(rs, "total_ms"),
                        "bytes": rs[0]["bytes"]}

    rs = sample(lambda: fetch(base, stream_path), a.runs)
    g.le("stream status!=200 count", sum(r["status"] != 200 for r in rs), 0, "")
    g.ge("stream prefixes", min(r.get("prefixes", 0) for r in rs), 1, "")
    g.le("stream TTFB median", med(rs, "ttfb_ms"), a.max_ttfb_ms)
    g.le("stream first-prefix median", med(rs, "first_prefix_ms") if all("first_prefix_ms" in r for r in rs) else float("inf"),
         a.max_first_prefix_ms)
    g.le("stream total median", med(rs, "total_ms"), a.max_stream_ms)
    g.le("stream total worst", worst(rs, "total_ms"), a.max_stream_ms * 2)
    g.ge("stream throughput median", med(rs, "mb_per_s"), a.min_stream_mbps, "MB/s")
    report["stream"] = {k: med(rs, k) for k in ("ttfb_ms", "first_prefix_ms", "total_ms", "mb_per_s", "prefixes_per_s", "bytes")
                        if all(k in r for r in rs)}

    gz = sample(lambda: fetch(base, stream_path, gzip_ok=True), a.runs)
    g.le("stream(gzip) status!=200 count", sum(r["status"] != 200 for r in gz), 0, "")
    g.le("stream(gzip) TTFB median", med(gz, "ttfb_ms"), a.max_ttfb_ms)
    g.le("stream(gzip) total median", med(gz, "total_ms"), a.max_stream_ms)
    report["stream_gzip"] = {k: med(gz, k) for k in ("ttfb_ms", "total_ms", "bytes")}

    cc = concurrent(base, stream_path, a.workers, a.rounds)
    g.le("concurrent non-200 responses", 0 if cc["all_ok"] else 1, 0, "")
    g.ge(f"concurrent x{a.workers} requests/s", cc["requests_per_s"], a.min_concurrent_rps, "req/s")
    g.le(f"concurrent x{a.workers} slowest", cc["max_total_ms"], a.max_concurrent_slowest_ms)
    report["concurrent"] = {"workers": a.workers, **cc}

    report["gates"] = [{"name": n, "value": v, "limit": l, "ok": ok} for n, v, l, _, ok in g.rows]
    report["ok"] = g.ok

    lines = [f"## Web viewer latency / throughput — trace {a.trace_id} (ε={a.epsilon}, δ={a.delta})", "",
             "| Metric | Value | Gate | |", "|---|---|---|---|"]
    for n, v, l, _, ok in g.rows:
        lines.append(f"| {n} | {v:.1f} | {l} | {'✅' if ok else '❌'} |")
    md = "\n".join(lines) + "\n"
    print(md)
    if a.json_out:
        with open(a.json_out, "w") as f:
            json.dump(report, f, indent=2)
    if a.markdown_out:
        with open(a.markdown_out, "w") as f:
            f.write(md)
    if not g.ok:
        print("❌ web perf gate failed", file=sys.stderr)
        return 1
    print("✅ web perf gates passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
