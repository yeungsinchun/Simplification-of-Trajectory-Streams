#!/usr/bin/env python3
"""Preliminary fair comparison: simplify vs DOTS, matched on Frechet error."""
import os, statistics, subprocess, re

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SIMPLIFY = os.path.join(REPO, "build-refactor", "simplify")
DOTS = os.path.join(REPO, "traj-compression", "online", "DOTS", "DOTS_adapted")
FRECHET = os.path.join(REPO, "scripts", "frechet")

REPEATS = 5
EPSILON = 0.5
SIMPLIFY_DELTAS = [5, 10, 20, 40, 80, 160]
DOTS_THRESHOLDS = [500, 2000, 10000, 50000, 120000, 300000]

def median_core_ms(fn):
    return statistics.median(fn() for _ in range(REPEATS))

def run_simplify(dsid, delta):
    def _run():
        res = subprocess.run([SIMPLIFY, str(dsid), "-e", str(EPSILON), "-d", str(delta)],
                             capture_output=True, text=True)
        m = re.search(r'^total\s+([0-9.]+)', res.stderr, re.MULTILINE)
        return float(m.group(1)) if m else 0.0
    core_ms = median_core_ms(_run)
    subprocess.run([SIMPLIFY, str(dsid), "-e", str(EPSILON), "-d", str(delta)],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    simp_path = os.path.join(REPO, "data", str(dsid), "simplify.txt")
    with open(simp_path) as f:
        n_pts = int(f.readline().strip())
    frechet_out = subprocess.check_output([FRECHET, "--id", str(dsid), "--batch", simp_path], text=True)
    m = re.search(r'simplify\.txt:\s*([0-9.]+)', frechet_out)
    frechet = float(m.group(1)) if m else 0.0
    return core_ms, n_pts, frechet

def run_dots(dsid, thr):
    def _run():
        res = subprocess.run([DOTS, str(dsid), str(thr)], capture_output=True, text=True)
        m = re.search(r'DOTS_CORE_MS\s+([0-9.]+)', res.stderr)
        return float(m.group(1)) if m else 0.0
    core_ms = median_core_ms(_run)
    subprocess.run([DOTS, str(dsid), str(thr)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    dots_path = os.path.join(REPO, "data", str(dsid), "dots_simplified.txt")
    with open(dots_path) as f:
        n_pts = int(f.readline().strip())
    frechet_out = subprocess.check_output([FRECHET, "--id", str(dsid), "--batch", dots_path], text=True)
    m = re.search(r'dots_simplified\.txt:\s*([0-9.]+)', frechet_out)
    frechet = float(m.group(1)) if m else 0.0
    return core_ms, n_pts, frechet

def compare(dsid):
    print(f"Dataset {dsid}")
    simp_runs = []
    for d in SIMPLIFY_DELTAS:
        try:
            core, pts, frech = run_simplify(dsid, d)
            simp_runs.append((d, core, pts, frech))
            print(f"  simplify d={d:3d}: core={core:7.2f}ms  pts={pts:4d}  frechet={frech:7.2f}")
        except Exception as e:
            print(f"  simplify d={d:3d}: FAIL ({e})")
    dots_runs = []
    for t in DOTS_THRESHOLDS:
        try:
            core, pts, frech = run_dots(dsid, t)
            dots_runs.append((t, core, pts, frech))
            print(f"  DOTS thr={t:6d}: core={core:7.2f}ms  pts={pts:4d}  frechet={frech:7.2f}")
        except Exception as e:
            print(f"  DOTS thr={t:6d}: FAIL ({e})")
    if not simp_runs or not dots_runs:
        print("  No valid runs to compare.\n")
        return
    simp_runs.sort(key=lambda x: x[3])
    dots_runs.sort(key=lambda x: x[3])
    lo = max(simp_runs[0][3], dots_runs[0][3])
    hi = min(simp_runs[-1][3], dots_runs[-1][3])
    if lo > hi:
        print("  No Frechet overlap, cannot match.\n")
        return
    target_f = (lo + hi) / 2.0
    s_pick = min(simp_runs, key=lambda x: abs(x[3] - target_f))
    d_pick = min(dots_runs, key=lambda x: abs(x[3] - target_f))
    speedup = s_pick[1] / d_pick[1] if d_pick[1] > 0 else float('inf')
    pts_ratio = d_pick[2] / s_pick[2] if s_pick[2] > 0 else float('inf')
    print(f"  MATCHED (target Frechet ~{target_f:.1f}, overlap [{lo:.1f}, {hi:.1f}]):")
    print(f"    simplify: {s_pick[1]:7.2f}ms  {s_pick[2]:4d} pts  Frechet={s_pick[3]:.2f}")
    print(f"    DOTS:     {d_pick[1]:7.2f}ms  {d_pick[2]:4d} pts  Frechet={d_pick[3]:.2f}")
    print(f"    -> DOTS core is {speedup:.0f}x faster; DOTS uses {pts_ratio:.2f}x the points for similar error\n")

if __name__ == "__main__":
    for ds in [50, 59, 61, 66, 92, 95, 103, 104]:
        compare(ds)
