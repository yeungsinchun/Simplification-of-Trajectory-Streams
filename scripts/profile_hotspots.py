#!/usr/bin/env python3
"""Profile simplify hotspots across epsilon/delta on taxi trajectories.

Runs build-release/simplify with --time, parses the hierarchical timing
summary plus counter lines, and writes:
  results/profile_sweep.csv   - one row per (id, epsilon, delta) with component ms
  results/profile_summary.txt - aggregated human-readable report
"""
from __future__ import annotations

import argparse
import csv
import re
import statistics
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / "build-release" / "simplify"
DATA = REPO / "data"
OUT_CSV = REPO / "results" / "profile_sweep.csv"
OUT_TXT = REPO / "results" / "profile_summary.txt"

CORE_RE = re.compile(r"SIMPLIFY_CORE_MS:\s*([\d.]+)")
COUNTER_RE = re.compile(
    r"^(BOUNDARY_CANDIDATES_SUM|STAB_STEPS|ALIVE_CANDIDATE_ITERS|"
    r"WEDGE_PRUNE_HITS|SIMPLIFIED_POINTS):\s*(\d+)"
)
# Matches a timing-summary data line: name, wall, self, calls, %total, %parent
TIMING_RE = re.compile(
    r"^(\s*)(\S.*?)\s+([\d.]+)\s+([\d.]+)\s+(\d+)\s+([\d.]+)%\s+([\d.]+)%\s*$"
)

HOTSPOTS = [
    "intersect",
    "find_F",
    "wedge_gi_disjoint",
    "get_longest_stab",
    "get_boundary_points",
    "get_conv_from_grid",
    "update_S",
]

EPSILONS = [0.25, 0.5, 1.0, 2.0, 5.0, 10.0]
DELTAS = [50, 100, 200, 400]


def count_points(path: Path) -> int:
    with path.open() as f:
        return int(f.readline().strip())


def parse_profile(stderr: str) -> dict:
    out: dict = {"core_ms": None, "hotspots": {}, "counters": {}}
    for line in stderr.splitlines():
        m = CORE_RE.match(line)
        if m:
            out["core_ms"] = float(m.group(1))
            continue
        m = COUNTER_RE.match(line)
        if m:
            out["counters"][m.group(1)] = int(m.group(2))
            continue
        m = TIMING_RE.match(line)
        if m:
            name = m.group(2).strip()
            wall = float(m.group(3))
            self_ms = float(m.group(4))
            calls = int(m.group(5))
            pct = float(m.group(6))
            out["hotspots"][name] = {
                "wall_ms": wall,
                "self_ms": self_ms,
                "calls": calls,
                "pct_total": pct,
            }
    return out


def run_once(idv: int, eps: float, delta: float, timeout: int) -> dict:
    cmd = [str(BIN), "--in", str(idv), "-e", str(eps), "-d", str(delta), "--time"]
    try:
        proc = subprocess.run(
            cmd, cwd=REPO, capture_output=True, text=True, timeout=timeout
        )
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "core_ms": None, "hotspots": {}, "counters": {}}
    parsed = parse_profile(proc.stderr)
    parsed["status"] = "ok" if proc.returncode == 0 and parsed["core_ms"] is not None else "error"
    parsed["stderr_tail"] = proc.stderr[-500:] if proc.returncode != 0 else ""
    return parsed


def median_run(idv: int, eps: float, delta: float, repeats: int, timeout: int) -> dict:
    runs = []
    for _ in range(repeats):
        r = run_once(idv, eps, delta, timeout)
        if r["status"] == "ok":
            runs.append(r)
    if not runs:
        return run_once(idv, eps, delta, timeout)
    # Pick the run whose core_ms is the median (so hotspot breakdown matches).
    runs.sort(key=lambda r: r["core_ms"])
    return runs[len(runs) // 2]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", type=int, default=1)
    ap.add_argument("--b", type=int, default=10)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--timeout", type=int, default=120)
    ap.add_argument("--epsilons", type=float, nargs="+", default=EPSILONS)
    ap.add_argument("--deltas", type=float, nargs="+", default=DELTAS)
    args = ap.parse_args()

    if not BIN.exists():
        print(f"missing binary: {BIN}", file=sys.stderr)
        return 1

    ids = [i for i in range(args.a, args.b + 1) if (DATA / str(i) / "original.txt").exists()]
    rows = []
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    print(f"Profiling {len(ids)} IDs x {len(args.epsilons)} eps x {len(args.deltas)} deltas "
          f"x {args.repeats} repeats")

    for idv in ids:
        npts = count_points(DATA / str(idv) / "original.txt")
        for eps in args.epsilons:
            for delta in args.deltas:
                r = median_run(idv, eps, delta, args.repeats, args.timeout)
                hs = r.get("hotspots", {})
                ctr = r.get("counters", {})
                row = {
                    "id": idv,
                    "original_points": npts,
                    "epsilon": eps,
                    "delta": delta,
                    "status": r["status"],
                    "core_ms": r.get("core_ms") if r.get("core_ms") is not None else -1,
                    "simplified_points": ctr.get("SIMPLIFIED_POINTS", -1),
                    "stab_steps": ctr.get("STAB_STEPS", -1),
                    "alive_iters": ctr.get("ALIVE_CANDIDATE_ITERS", -1),
                    "wedge_prune_hits": ctr.get("WEDGE_PRUNE_HITS", -1),
                    "boundary_candidates_sum": ctr.get("BOUNDARY_CANDIDATES_SUM", -1),
                }
                for name in HOTSPOTS:
                    h = hs.get(name, {})
                    row[f"{name}_ms"] = h.get("wall_ms", -1)
                    row[f"{name}_pct"] = h.get("pct_total", -1)
                    row[f"{name}_calls"] = h.get("calls", -1)
                # get_longest_stab self time is the unaccounted loop overhead
                gls = hs.get("get_longest_stab", {})
                row["get_longest_stab_self_ms"] = gls.get("self_ms", -1)
                rows.append(row)
                print(
                    f"[id={idv:2d} n={npts:5d}] e={eps:5.2f} d={delta:6.1f}  "
                    f"core={row['core_ms']:8.2f} ms  "
                    f"intersect={row.get('intersect_pct', -1):5.1f}%  "
                    f"find_F={row.get('find_F_pct', -1):5.1f}%  "
                    f"wedge={row.get('wedge_gi_disjoint_pct', -1):5.1f}%  "
                    f"simp={row['simplified_points']:4d}  [{row['status']}]"
                )

    fieldnames = list(rows[0].keys()) if rows else []
    with OUT_CSV.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    # ---- aggregate report ----
    ok = [r for r in rows if r["status"] == "ok"]
    lines = []
    lines.append("PROFILE SUMMARY - streaming delta-simplification on T-Drive taxi data")
    lines.append("=" * 72)
    lines.append(f"IDs: {ids}")
    lines.append(f"epsilons: {args.epsilons}")
    lines.append(f"deltas: {args.deltas}")
    lines.append(f"repeats/median: {args.repeats}")
    lines.append(f"ok rows: {len(ok)} / {len(rows)}")
    lines.append("")

    # Hotspot share averaged over all ok runs
    lines.append("1) Where does time go? (mean % of core across all ok runs)")
    lines.append("-" * 72)
    for name in HOTSPOTS:
        vals = [r[f"{name}_pct"] for r in ok if r[f"{name}_pct"] >= 0]
        if vals:
            lines.append(f"  {name:24s}  mean={statistics.mean(vals):5.1f}%  "
                         f"median={statistics.median(vals):5.1f}%  "
                         f"max={max(vals):5.1f}%")
    self_vals = [r["get_longest_stab_self_ms"] / r["core_ms"] * 100
                 for r in ok if r["core_ms"] > 0 and r["get_longest_stab_self_ms"] >= 0]
    if self_vals:
        lines.append(f"  {'get_longest_stab self':24s}  mean={statistics.mean(self_vals):5.1f}%  "
                     f"median={statistics.median(self_vals):5.1f}%  "
                     f"max={max(self_vals):5.1f}%")
    lines.append("")

    # Epsilon effect at fixed delta=200 (or nearest available)
    fixed_d = 200.0 if 200.0 in args.deltas else args.deltas[len(args.deltas) // 2]
    lines.append(f"2) Epsilon sweep at delta={fixed_d} (median core_ms across IDs)")
    lines.append("-" * 72)
    lines.append(f"  {'eps':>6}  {'core_ms':>10}  {'intersect%':>10}  {'find_F%':>8}  "
                 f"{'alive_iters':>12}  {'simp_pts':>8}")
    for eps in args.epsilons:
        subset = [r for r in ok if r["epsilon"] == eps and r["delta"] == fixed_d]
        if not subset:
            continue
        lines.append(
            f"  {eps:6.2f}  {statistics.median(r['core_ms'] for r in subset):10.2f}  "
            f"{statistics.median(r['intersect_pct'] for r in subset):10.1f}  "
            f"{statistics.median(r['find_F_pct'] for r in subset):8.1f}  "
            f"{statistics.median(r['alive_iters'] for r in subset):12.0f}  "
            f"{statistics.median(r['simplified_points'] for r in subset):8.0f}"
        )
    lines.append("")

    # Delta effect at fixed epsilon=0.5
    fixed_e = 0.5 if 0.5 in args.epsilons else args.epsilons[0]
    lines.append(f"3) Delta sweep at epsilon={fixed_e} (median core_ms across IDs)")
    lines.append("-" * 72)
    lines.append(f"  {'delta':>6}  {'core_ms':>10}  {'intersect%':>10}  {'find_F%':>8}  "
                 f"{'stab_steps':>10}  {'simp_pts':>8}")
    for delta in args.deltas:
        subset = [r for r in ok if r["epsilon"] == fixed_e and r["delta"] == delta]
        if not subset:
            continue
        lines.append(
            f"  {delta:6.1f}  {statistics.median(r['core_ms'] for r in subset):10.2f}  "
            f"{statistics.median(r['intersect_pct'] for r in subset):10.1f}  "
            f"{statistics.median(r['find_F_pct'] for r in subset):8.1f}  "
            f"{statistics.median(r['stab_steps'] for r in subset):10.0f}  "
            f"{statistics.median(r['simplified_points'] for r in subset):8.0f}"
        )
    lines.append("")

    # Per-ID at default-ish params
    lines.append(f"4) Per-trajectory at e={fixed_e}, d={fixed_d}")
    lines.append("-" * 72)
    lines.append(f"  {'id':>4}  {'n':>6}  {'core_ms':>10}  {'intersect%':>10}  "
                 f"{'find_F%':>8}  {'alive':>10}  {'simp':>6}")
    for idv in ids:
        subset = [r for r in ok if r["id"] == idv and r["epsilon"] == fixed_e and r["delta"] == fixed_d]
        if not subset:
            continue
        r = subset[0]
        lines.append(
            f"  {idv:4d}  {r['original_points']:6d}  {r['core_ms']:10.2f}  "
            f"{r['intersect_pct']:10.1f}  {r['find_F_pct']:8.1f}  "
            f"{r['alive_iters']:10d}  {r['simplified_points']:6d}"
        )
    lines.append("")

    # Scaling: core_ms vs alive_iters correlation hint
    lines.append("5) Cost drivers")
    lines.append("-" * 72)
    if ok:
        by_alive = sorted(ok, key=lambda r: r["alive_iters"])
        lines.append(
            f"  alive_candidate_iters range: "
            f"{by_alive[0]['alive_iters']} .. {by_alive[-1]['alive_iters']}"
        )
        lines.append(
            f"  core_ms range: "
            f"{min(r['core_ms'] for r in ok):.2f} .. {max(r['core_ms'] for r in ok):.2f}"
        )
        # Pn proxy: boundary_candidates_sum / number of prefixes ~= avg Pn
        # prefixes ~= simplified_points/2
        pn_est = []
        for r in ok:
            prefixes = max(1, r["simplified_points"] // 2)
            if r["boundary_candidates_sum"] > 0:
                pn_est.append(r["boundary_candidates_sum"] / prefixes)
        if pn_est:
            lines.append(
                f"  estimated avg boundary candidates |P|: "
                f"median={statistics.median(pn_est):.0f}  "
                f"max={max(pn_est):.0f}"
            )
        prune_rates = [
            r["wedge_prune_hits"] / r["alive_iters"]
            for r in ok if r["alive_iters"] > 0 and r["wedge_prune_hits"] >= 0
        ]
        if prune_rates:
            lines.append(
                f"  wedge prune hit rate: "
                f"median={100 * statistics.median(prune_rates):.1f}%  "
                f"max={100 * max(prune_rates):.1f}%"
            )
    lines.append("")
    lines.append(f"Raw CSV: {OUT_CSV}")

    report = "\n".join(lines) + "\n"
    OUT_TXT.write_text(report)
    print()
    print(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
