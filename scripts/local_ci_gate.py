#!/usr/bin/env python3
"""Local CI gate: compare NEW simplify binary vs BASELINE on IDs 1..10.

Mirrors .github/workflows/correctness.yml and benchmark.yml:
  - Frechet: new_dist <= orig_dist + DIST_TOL
  - Points:  new_points <= orig_points + POINTS_TOL
  - Perf:    new_ms <= orig_ms * SLOWDOWN_LIMIT  (when orig_ms >= MIN_BENCH_MS)

Uses the same CI params (EPSILON=299, DELTA=1) by default.
Also supports a 'profile' mode with e=0.5 d=200 for hotspot-relevant checks.
"""
from __future__ import annotations

import argparse
import re
import statistics
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CORE_RE = re.compile(r"SIMPLIFY_CORE_MS:\s*([\d.]+)")
RAW_DIST_RE = re.compile(r"^[0-9]+(?:\.[0-9]+)?$", re.M)
LABEL_DIST_RE = re.compile(r"Frechet distance:\s*([0-9]+(?:\.[0-9]+)?)", re.I)


def count_points(path: Path) -> int:
    with path.open() as f:
        return int(f.readline().strip())


def run_simplify(bin_path: Path, idv: int, eps: float, delta: float,
                 *, dist: bool, cwd: Path, timeout: int = 180) -> tuple[str, str]:
    cmd = [str(bin_path), str(idv), "-e", str(eps), "-d", str(delta)]
    if dist:
        cmd.append("--dist")
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
    if proc.returncode != 0:
        raise RuntimeError(
            f"{bin_path} id={idv} failed rc={proc.returncode}\n"
            f"stdout:\n{proc.stdout[-1000:]}\nstderr:\n{proc.stderr[-1000:]}"
        )
    return proc.stdout, proc.stderr


def parse_core_ms(stderr: str, stdout: str) -> float:
    for text in (stderr, stdout):
        m = CORE_RE.search(text)
        if m:
            return float(m.group(1))
    raise RuntimeError("missing SIMPLIFY_CORE_MS")


def parse_dist(stdout: str, stderr: str) -> float:
    text = stdout + "\n" + stderr
    ms = RAW_DIST_RE.findall(text)
    if ms:
        return float(ms[-1])
    m = LABEL_DIST_RE.search(text)
    if m:
        return float(m.group(1))
    raise RuntimeError(f"could not parse Frechet distance\n{text[-1500:]}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", type=Path, required=True)
    ap.add_argument("--new", type=Path, required=True)
    ap.add_argument("--baseline-cwd", type=Path, default=None,
                    help="cwd for baseline runs (isolates data/<id>/simplify.txt)")
    ap.add_argument("--new-cwd", type=Path, default=None)
    ap.add_argument("--a", type=int, default=1)
    ap.add_argument("--b", type=int, default=10)
    ap.add_argument("--epsilon", type=float, default=299.0)
    ap.add_argument("--delta", type=float, default=1.0)
    ap.add_argument("--dist-tol", type=float, default=0.01)
    ap.add_argument("--points-tol", type=int, default=0)
    ap.add_argument("--slowdown-limit", type=float, default=1.50)
    ap.add_argument("--min-bench-ms", type=float, default=20.0)
    ap.add_argument("--bench-runs", type=int, default=3)
    ap.add_argument("--skip-dist", action="store_true")
    ap.add_argument("--skip-bench", action="store_true")
    args = ap.parse_args()

    base_cwd = args.baseline_cwd or REPO
    new_cwd = args.new_cwd or REPO
    ids = list(range(args.a, args.b + 1))

    overall_ok = True
    print(f"Gate params: e={args.epsilon} d={args.delta} "
          f"dist_tol={args.dist_tol} points_tol={args.points_tol} "
          f"slowdown={args.slowdown_limit}x min_ms={args.min_bench_ms}")
    print(f"baseline={args.baseline}")
    print(f"new     ={args.new}")
    print()

    # ---- Correctness ----
    if not args.skip_dist:
        print("=== Correctness (Frechet + points) ===")
        print(f"{'id':>4}  {'orig_dist':>12}  {'new_dist':>12}  {'d_dist':>10}  "
              f"{'orig_pts':>8}  {'new_pts':>8}  {'d_pts':>6}  status")
        for idv in ids:
            out_b, err_b = run_simplify(args.baseline, idv, args.epsilon, args.delta,
                                        dist=True, cwd=base_cwd)
            orig_dist = parse_dist(out_b, err_b)
            orig_pts = count_points(base_cwd / "data" / str(idv) / "simplify.txt")

            out_n, err_n = run_simplify(args.new, idv, args.epsilon, args.delta,
                                        dist=True, cwd=new_cwd)
            new_dist = parse_dist(out_n, err_n)
            new_pts = count_points(new_cwd / "data" / str(idv) / "simplify.txt")

            status = "OK"
            if new_dist > orig_dist + args.dist_tol + 1e-9:
                status = "FAIL_DIST"
                overall_ok = False
            elif new_pts - orig_pts > args.points_tol:
                status = "FAIL_POINTS"
                overall_ok = False

            print(f"{idv:4d}  {orig_dist:12.6f}  {new_dist:12.6f}  "
                  f"{new_dist - orig_dist:10.6f}  {orig_pts:8d}  {new_pts:8d}  "
                  f"{new_pts - orig_pts:6d}  {status}")
        print()

    # ---- Benchmark ----
    if not args.skip_bench:
        print(f"=== Benchmark (avg of {args.bench_runs}, limit {args.slowdown_limit}x) ===")
        print(f"{'id':>4}  {'prev_ms':>10}  {'new_ms':>10}  {'ratio':>8}  status")
        for idv in ids:
            orig_times = []
            new_times = []
            for _ in range(args.bench_runs):
                _, err = run_simplify(args.baseline, idv, args.epsilon, args.delta,
                                      dist=False, cwd=base_cwd)
                orig_times.append(parse_core_ms(err, ""))
                _, err = run_simplify(args.new, idv, args.epsilon, args.delta,
                                      dist=False, cwd=new_cwd)
                new_times.append(parse_core_ms(err, ""))
            orig_ms = statistics.mean(orig_times)
            new_ms = statistics.mean(new_times)
            ratio = new_ms / orig_ms if orig_ms > 0 else float("inf")
            status = "OK"
            if orig_ms >= args.min_bench_ms and new_ms > orig_ms * args.slowdown_limit + 1e-9:
                status = "FAIL"
                overall_ok = False
            print(f"{idv:4d}  {orig_ms:10.4f}  {new_ms:10.4f}  {ratio:8.3f}  {status}")
        print()

    if overall_ok:
        print("✅ Local CI gate OK")
        return 0
    print("❌ Local CI gate FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
