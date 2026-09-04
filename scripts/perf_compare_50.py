#!/usr/bin/env python3
"""Compare baseline vs optimized simplify core_ms over many taxi IDs.

Fixed (epsilon, delta) pairs. Writes CSV + PNG plots under results/.
"""
from __future__ import annotations

import argparse
import csv
import re
import statistics
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt

CORE_RE = re.compile(r"SIMPLIFY_CORE_MS:\s*([\d.]+)")


def count_points(path: Path) -> int:
    with path.open() as f:
        return int(f.readline().strip())


def core_ms(bin_path: Path, cwd: Path, idv: int, eps: float, delta: float,
            repeats: int) -> float:
    times = []
    for _ in range(repeats):
        proc = subprocess.run(
            [str(bin_path), str(idv), "-e", str(eps), "-d", str(delta)],
            cwd=cwd, capture_output=True, text=True, timeout=300,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"{bin_path} id={idv} e={eps} d={delta} rc={proc.returncode}\n"
                f"{proc.stderr[-500:]}"
            )
        m = CORE_RE.search(proc.stderr) or CORE_RE.search(proc.stdout)
        if not m:
            raise RuntimeError(f"missing SIMPLIFY_CORE_MS for id={idv}")
        times.append(float(m.group(1)))
    return statistics.median(times)


def summarize(ratios: list[float]) -> dict:
    return {
        "n": len(ratios),
        "mean": statistics.mean(ratios),
        "median": statistics.median(ratios),
        "stdev": statistics.stdev(ratios) if len(ratios) > 1 else 0.0,
        "min": min(ratios),
        "max": max(ratios),
        "pct_faster": 100.0 * sum(1 for r in ratios if r < 0.98) / len(ratios),
        "pct_slower": 100.0 * sum(1 for r in ratios if r > 1.02) / len(ratios),
        "pct_flat": 100.0 * sum(1 for r in ratios if 0.98 <= r <= 1.02) / len(ratios),
        "geomean": statistics.geometric_mean(ratios),
    }


def plot_pair(rows: list[dict], eps: float, delta: float, out_dir: Path) -> Path:
    rows = sorted(rows, key=lambda r: r["id"])
    ids = [r["id"] for r in rows]
    base = [r["base_ms"] for r in rows]
    opt = [r["opt_ms"] for r in rows]
    ratios = [r["ratio"] for r in rows]
    npts = [r["n_points"] for r in rows]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(f"Baseline vs optimized simplify (e={eps}, d={delta}, N={len(rows)} IDs)",
                 fontsize=13)

    ax = axes[0, 0]
    ax.plot(ids, base, "o-", label="baseline (main)", markersize=3, linewidth=1)
    ax.plot(ids, opt, "s-", label="optimized (P0+P1)", markersize=3, linewidth=1)
    ax.set_xlabel("trajectory id")
    ax.set_ylabel("core_ms (median)")
    ax.set_title("Wall time by ID")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    colors = ["#2ca02c" if r < 1 else "#d62728" for r in ratios]
    ax.bar(ids, ratios, color=colors, width=0.8)
    ax.axhline(1.0, color="black", linewidth=1)
    ax.axhline(1.5, color="gray", linestyle="--", linewidth=0.8, label="CI 1.5x bar")
    ax.set_xlabel("trajectory id")
    ax.set_ylabel("opt / baseline")
    ax.set_title("Speedup ratio (<1 is faster)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis="y")

    ax = axes[1, 0]
    ax.scatter(npts, ratios, c=colors, s=28, alpha=0.85)
    ax.axhline(1.0, color="black", linewidth=1)
    ax.set_xlabel("original points")
    ax.set_ylabel("opt / baseline")
    ax.set_title("Ratio vs trajectory length")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.hist(ratios, bins=min(15, max(5, len(ratios) // 3)), color="#1f77b4", edgecolor="white")
    ax.axvline(1.0, color="black", linewidth=1)
    ax.axvline(statistics.median(ratios), color="#ff7f0e", linestyle="--",
               label=f"median={statistics.median(ratios):.3f}")
    ax.axvline(statistics.mean(ratios), color="#2ca02c", linestyle=":",
               label=f"mean={statistics.mean(ratios):.3f}")
    ax.set_xlabel("opt / baseline")
    ax.set_ylabel("count")
    ax.set_title("Ratio distribution")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis="y")

    fig.tight_layout()
    tag = f"e{eps}_d{delta}".replace(".", "p")
    path = out_dir / f"perf_compare_50_{tag}.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline-bin", type=Path, required=True)
    ap.add_argument("--baseline-cwd", type=Path, required=True)
    ap.add_argument("--opt-bin", type=Path, required=True)
    ap.add_argument("--opt-cwd", type=Path, required=True)
    ap.add_argument("--a", type=int, default=1)
    ap.add_argument("--b", type=int, default=50)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--epsilon", type=float, nargs="+", default=[0.5, 299.0])
    ap.add_argument("--delta", type=float, nargs="+", default=[200.0, 1.0])
    ap.add_argument("--out-dir", type=Path, default=Path("results"))
    args = ap.parse_args()

    if len(args.epsilon) != len(args.delta):
        print("epsilon and delta lists must have the same length", file=sys.stderr)
        return 1

    args.out_dir.mkdir(parents=True, exist_ok=True)
    ids = list(range(args.a, args.b + 1))

    all_summaries = []
    for eps, delta in zip(args.epsilon, args.delta):
        print(f"\n=== e={eps} d={delta} over IDs {args.a}..{args.b} "
              f"({args.repeats} repeats, median) ===", flush=True)
        rows = []
        for idv in ids:
            orig = args.opt_cwd / "data" / str(idv) / "original.txt"
            if not orig.exists():
                print(f"  skip id={idv}: missing original.txt", flush=True)
                continue
            npts = count_points(orig)
            b = core_ms(args.baseline_bin, args.baseline_cwd, idv, eps, delta, args.repeats)
            o = core_ms(args.opt_bin, args.opt_cwd, idv, eps, delta, args.repeats)
            ratio = o / b if b > 0 else float("inf")
            row = {
                "id": idv,
                "n_points": npts,
                "epsilon": eps,
                "delta": delta,
                "base_ms": round(b, 4),
                "opt_ms": round(o, 4),
                "ratio": round(ratio, 4),
                "saved_ms": round(b - o, 4),
            }
            rows.append(row)
            print(
                f"  id={idv:3d} n={npts:5d}  base={b:8.3f}  opt={o:8.3f}  "
                f"ratio={ratio:6.3f}",
                flush=True,
            )

        tag = f"e{eps}_d{delta}".replace(".", "p")
        csv_path = args.out_dir / f"perf_compare_50_{tag}.csv"
        with csv_path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

        ratios = [r["ratio"] for r in rows]
        stats = summarize(ratios)
        stats["epsilon"] = eps
        stats["delta"] = delta
        stats["total_base_ms"] = sum(r["base_ms"] for r in rows)
        stats["total_opt_ms"] = sum(r["opt_ms"] for r in rows)
        stats["total_saved_ms"] = stats["total_base_ms"] - stats["total_opt_ms"]
        all_summaries.append(stats)

        png = plot_pair(rows, eps, delta, args.out_dir)
        print(f"CSV: {csv_path}")
        print(f"PNG: {png}")
        print(
            f"STATS e={eps} d={delta}: n={stats['n']}  "
            f"mean_ratio={stats['mean']:.3f}  median_ratio={stats['median']:.3f}  "
            f"geomean={stats['geomean']:.3f}  "
            f"faster={stats['pct_faster']:.0f}%  slower={stats['pct_slower']:.0f}%  "
            f"flat={stats['pct_flat']:.0f}%  "
            f"total_base={stats['total_base_ms']:.1f}ms  "
            f"total_opt={stats['total_opt_ms']:.1f}ms  "
            f"saved={stats['total_saved_ms']:.1f}ms"
        )

    summary_path = args.out_dir / "perf_compare_50_summary.txt"
    with summary_path.open("w") as f:
        f.write("Performance comparison: baseline (main) vs optimized (wedge prune + Gi dedup/skip to_ccw)\n")
        f.write(f"IDs {args.a}-{args.b}, repeats={args.repeats} (median)\n\n")
        for s in all_summaries:
            f.write(f"e={s['epsilon']} d={s['delta']}\n")
            f.write(f"  n={s['n']}\n")
            f.write(f"  mean_ratio={s['mean']:.4f}\n")
            f.write(f"  median_ratio={s['median']:.4f}\n")
            f.write(f"  geomean_ratio={s['geomean']:.4f}\n")
            f.write(f"  stdev={s['stdev']:.4f}  min={s['min']:.4f}  max={s['max']:.4f}\n")
            f.write(f"  pct_faster(<0.98)={s['pct_faster']:.1f}%\n")
            f.write(f"  pct_slower(>1.02)={s['pct_slower']:.1f}%\n")
            f.write(f"  pct_flat=[0.98,1.02]={s['pct_flat']:.1f}%\n")
            f.write(f"  total_base_ms={s['total_base_ms']:.2f}\n")
            f.write(f"  total_opt_ms={s['total_opt_ms']:.2f}\n")
            f.write(f"  total_saved_ms={s['total_saved_ms']:.2f}\n\n")
    print(f"\nSummary: {summary_path}")
    print(summary_path.read_text())
    return 0


if __name__ == "__main__":
    sys.exit(main())
