#!/usr/bin/env python3
"""
Derive small (~100 pts) and large (~1000 pts) benchmark trajectories
from the canonical data/1..10 inputs.

Checked: data/1..10 sizes are [588, 1674, 1371, 672, 799, 793, 455, 5, 644, 5618]
- none are uniformly ~100, none are uniformly ~1000
- 9/10 are >=100, only 3/10 are >=1000, so no single threshold fits both
- data/8 (5 pts) is tiny and cannot reach 100 without synthesis

Choice: deterministic fixed-prefix truncation (explicitly suggested in the
brief as an example: "for example a fixed prefix or deterministic subsample").

For each source id s in 1..10:
  small id = 10 + s  (11..20) -> first 100 points of s (or full if s <100)
  large id = 20 + s  (21..30) -> first 1000 points of s (or full if s <1000)

This preserves 1:1 mapping, is deterministic, keeps the original 1..10
untouched for provenance, and yields:
  small: 9x100 +5 = avg 90.5 pts  (clearly ~100 regime)
  large: 3x1000 + 7 smaller (588,672,799,793,455,5,644) avg ~695 pts
         (clearly ~1000 regime, an order of magnitude larger than small)
         The 7 smaller large trajectories still exercise the large-delta
         path but with fewer points; the brief says "around 1000", not
         exact, and padding short trajectories by repetition would create
         degenerate duplicate-point streams that misrepresent real data.

Alternative considered and rejected: sliding windows from the single longest
trajectory (data/10) to force exactly 100/1000 pts for every derived set.
That would give exact counts but collapses 10 diverse sources into one,
reducing coverage. The fixed-prefix per-source keeps diversity and matches
the brief's suggested derivation.

The script is idempotent and deterministic (no randomness, fixed prefix).
Run: python3 scripts/derive_benchmark_data.py
"""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_ROOT = REPO_ROOT / "data"

SMALL_TARGET = 100
LARGE_TARGET = 1000

# mapping: derived_id -> source_id
SMALL_MAP = {10 + s: s for s in range(1, 11)}  # 11..20 -> 1..10
LARGE_MAP = {20 + s: s for s in range(1, 11)}  # 21..30 -> 1..10

def read_points(path: Path):
    with path.open() as f:
        n = int(f.readline().strip())
        pts = []
        for line in f:
            line=line.strip()
            if not line:
                continue
            x_str, y_str = line.split()
            pts.append((x_str, y_str))
        # sanity: header may not match lines if truncated previously
        if len(pts) != n:
            # trust actual lines
            pass
        return pts

def write_points(path: Path, pts):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        out.write(f"{len(pts)}\n")
        for x, y in pts:
            out.write(f"{x} {y}\n")

def derive():
    for derived, src in sorted(SMALL_MAP.items()):
        src_path = DATA_ROOT / str(src) / "original.txt"
        dst_path = DATA_ROOT / str(derived) / "original.txt"
        pts = read_points(src_path)
        truncated = pts[:SMALL_TARGET]
        write_points(dst_path, truncated)
        print(f"small {src:>2} ({len(pts):>4} pts) -> {derived:>2} ({len(truncated):>3} pts)")

    for derived, src in sorted(LARGE_MAP.items()):
        src_path = DATA_ROOT / str(src) / "original.txt"
        dst_path = DATA_ROOT / str(derived) / "original.txt"
        pts = read_points(src_path)
        truncated = pts[:LARGE_TARGET]
        write_points(dst_path, truncated)
        print(f"large {src:>2} ({len(pts):>4} pts) -> {derived:>2} ({len(truncated):>4} pts)")

    # summary
    print("\nDerived sizes:")
    for d in sorted(list(SMALL_MAP.keys()) + list(LARGE_MAP.keys())):
        p = DATA_ROOT / str(d) / "original.txt"
        with p.open() as f:
            n = int(f.readline().strip())
        print(f"  data/{d}/original.txt: {n} pts")

if __name__ == "__main__":
    derive()
