#!/usr/bin/env python3
"""Summarize the benchmark CSV: timing and Frechet-error by baseline.

Usage:
    python3 summarize.py [CSV_PATH]

Default path is the repo-root compare_points.csv (which is where benchmark.py
writes its output). It computes per-baseline means/medians of:
  - baseline_time, simplify_time
  - baseline_frechet, simplify_frechet
  - simplify_points / baseline_points ratio (compression)
  - speedup (baseline_time vs simplify_time)
"""
import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path


def fmt_num(x, w=10, prec=2):
    if x is None:
        return " " * w
    return f"{x:>{w}.{prec}f}"


def main():
    csv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else (
        Path(__file__).resolve().parent.parent / "compare_points.csv"
    )
    if not csv_path.exists():
        print(f"CSV not found: {csv_path}", file=sys.stderr)
        sys.exit(1)

    by_baseline = defaultdict(list)
    with csv_path.open() as f:
        for row in csv.DictReader(f):
            try:
                algo = row["baseline_algo"]
                bp = float(row["baseline_points"])
                sp = float(row["simplify_points"])
                bf = float(row["baseline_frechet"])
                sf = float(row["simplify_frechet"])
                bt = float(row["baseline_time"])
                st = float(row["simplify_time"])
                if sp <= 0 or bp <= 0 or sf <= 0 or bf <= 0:
                    continue
                by_baseline[algo].append({
                    "id": int(row["id"]),
                    "ratio": sp / bp,
                    "fret_ratio": sf / bf,
                    "speedup": bt / st if st > 0 else float("inf"),
                    "simplify_time": st,
                    "baseline_time": bt,
                })
            except (KeyError, ValueError):
                continue

    print(f"Loaded {sum(len(v) for v in by_baseline.values())} rows from {csv_path}\n")
    header = f"{'baseline':<10} {'n':>4}  {'pts_ratio':>10}  {'fret_ratio':>11}  {'simp_time':>10}  {'base_time':>10}  {'speedup':>10}"
    print(header)
    print("-" * len(header))
    for algo, rows in sorted(by_baseline.items()):
        n = len(rows)
        med_ratio = statistics.median(r["ratio"] for r in rows)
        med_fret = statistics.median(r["fret_ratio"] for r in rows)
        med_st = statistics.median(r["simplify_time"] for r in rows)
        med_bt = statistics.median(r["baseline_time"] for r in rows)
        med_su = statistics.median(r["speedup"] for r in rows)
        print(
            f"{algo:<10} {n:>4}  {med_ratio:>10.3f}  {med_fret:>11.4f}  "
            f"{med_st:>10.3f}  {med_bt:>10.4f}  {med_su:>10.2f}"
        )
    print("\nAccuracy (fret_ratio <= 1.01 means simplify at least as accurate):")
    for algo, rows in sorted(by_baseline.items()):
        ok = sum(1 for r in rows if r["fret_ratio"] <= 1.01)
        print(f"  {algo:<10}: {ok}/{len(rows)} ({100*ok/max(1,len(rows)):.0f}%)")
    print("\nDefinitions:")
    print("  pts_ratio    = median(simplify_points / baseline_points)")
    print("                  (>1 means simplify keeps MORE points than the baseline)")
    print("  fret_ratio   = median(simplify_frechet / baseline_frechet)")
    print("                  (>1 means simplify is less accurate than the baseline)")
    print("  speedup      = median(baseline_time / simplify_time)")
    print("                  (>1 means simplify is FASTER than the baseline)")
    print("  simp_time    = median(simplify's algorithm-only wall time in seconds)")
    print("  base_time    = median(baseline binary's wall time in seconds)")


if __name__ == "__main__":
    main()