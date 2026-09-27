#!/usr/bin/env python3
"""Generate charts for SOTS vs SQUISH benchmark."""
import argparse
import json
import pathlib
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

REPO = pathlib.Path(__file__).resolve().parent.parent
OUT_DIR = REPO / ".lavish" / "charts"
OUT_DIR.mkdir(parents=True, exist_ok=True)

parser = argparse.ArgumentParser(description="Generate SOTS vs SQUISH charts")
parser.add_argument("--input", type=pathlib.Path, help="Path to bench-squish JSON (default: latest large)")
args = parser.parse_args()

if args.input:
    JSON_PATH = args.input
else:
    candidates = sorted((REPO / ".lavish").glob("bench-squish-*.json"))
    if not candidates:
        print("No bench-squish JSON found", file=sys.stderr)
        sys.exit(1)
    large_candidates = []
    for p in candidates:
        try:
            with open(p) as f:
                j = json.load(f)
            if j.get("meta", {}).get("size") == "large":
                large_candidates.append(p)
        except Exception:
            continue
    if large_candidates:
        JSON_PATH = large_candidates[-1]
    else:
        print("No large bench-squish JSON found; use --input to specify", file=sys.stderr)
        sys.exit(1)

print(f"Reading {JSON_PATH}")
with open(JSON_PATH) as f:
    data = json.load(f)

aggs = data["aggregates"]
order = {299: 0, 30: 1, 5: 2, 0.5: 3, 0.1: 4}
aggs_sorted = sorted(aggs, key=lambda x: order.get(x["epsilon"], 99))

labels = []
tier_names = {299: "extra-coarse\n299", 30: "coarse\n30", 5: "mid\n5", 0.5: "fine\n0.5", 0.1: "extra-fine\n0.1"}
for a in aggs_sorted:
    labels.append(tier_names.get(a["epsilon"], str(a["epsilon"])))

sots_means = [a["mean_sots_ms"] for a in aggs_sorted]
squish_means = [a["mean_squish_ms"] for a in aggs_sorted]

speedups = [a["mean_speedup_squish_over_sots"] for a in aggs_sorted]
overheads = [1/s if s else 0 for s in speedups]

fig, ax = plt.subplots(figsize=(10, 6))
x = range(len(labels))
width = 0.35
ax.bar([i - width/2 for i in x], sots_means, width, label="SOTS (simplify)", color="#2563eb")
ax.bar([i + width/2 for i in x], squish_means, width, label="SQUISH", color="#16a34a")
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel("Mean core time (ms, log scale)")
ax.set_yscale("log")
ax.set_title("SOTS vs SQUISH — Mean core time (large datasets, 10 runs)")
ax.legend()
ax.grid(True, which="both", alpha=0.3)
for i, (s, q) in enumerate(zip(sots_means, squish_means)):
    ax.text(i - width/2, s*1.1, f"{s:.1f}", ha="center", va="bottom", fontsize=8, color="#2563eb")
    ax.text(i + width/2, q*1.1, f"{q:.2f}", ha="center", va="bottom", fontsize=8, color="#16a34a")

plt.tight_layout()
out1 = OUT_DIR / "runtime_bars.png"
plt.savefig(out1, dpi=180)
print(f"Wrote {out1}")
plt.close()

fig, ax = plt.subplots(figsize=(10, 6))
bars = ax.bar(x, speedups, width=0.6, color="#9333ea", edgecolor="black")
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel("Speedup (SQUISH / SOTS) — <1 means SQUISH faster")
ax.set_title("Speedup vs SQUISH (large, 10-run means) — SOTS is slower, SQUISH wins")
ax.set_ylim(0, 0.06)
ax.grid(True, axis="y", alpha=0.3)
for i, v in enumerate(speedups):
    ax.text(i, v + 0.001, f"{v:.3f}×", ha="center", va="bottom", fontsize=9, fontweight="bold")
    overhead = overheads[i]
    ax.text(i, 0.005, f"({overhead:.0f}× slower)", ha="center", va="bottom", fontsize=7, color="#555")

plt.tight_layout()
out2 = OUT_DIR / "speedup_bars.png"
plt.savefig(out2, dpi=180)
print(f"Wrote {out2}")
plt.close()

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(labels, sots_means, marker="o", linewidth=3, color="#2563eb", label="SOTS")
ax.plot(labels, squish_means, marker="s", linewidth=3, color="#16a34a", label="SQUISH")
ax.set_yscale("log")
ax.set_ylabel("Mean core time (ms, log scale)")
ax.set_title("Runtime vs epsilon tier (large, delta=300/(1+ε))")
ax.legend()
ax.grid(True, which="both", alpha=0.3)
ax.annotate(f"SOTS {sots_means[-1]:.0f}ms vs SQUISH {squish_means[-1]:.2f}ms\n({overheads[-1]:.0f}×)", xy=(4, sots_means[-1]), xytext=(2.5, 50),
            arrowprops=dict(arrowstyle="->", color="black"), fontsize=9, ha="center",
            bbox=dict(boxstyle="round,pad=0.3", fc="yellow", alpha=0.5))

plt.tight_layout()
out3 = OUT_DIR / "runtime_lines.png"
plt.savefig(out3, dpi=180)
print(f"Wrote {out3}")
plt.close()

import collections
points_by_eps = collections.defaultdict(list)
squish_points_by_eps = collections.defaultdict(list)
for c in data["cases"]:
    points_by_eps[c["epsilon"]].append(c["sots_points"])
    squish_points_by_eps[c["epsilon"]].append(c["squish_points"])

fig, ax = plt.subplots(figsize=(10, 5))
sots_pts_means = [sum(points_by_eps[e])/len(points_by_eps[e]) for e in [299,30,5,0.5,0.1]]
squish_pts_means = [sum(squish_points_by_eps[e])/len(squish_points_by_eps[e]) for e in [299,30,5,0.5,0.1]]
ax.bar([i - width/2 for i in x], sots_pts_means, width, label="SOTS points kept", color="#2563eb")
ax.bar([i + width/2 for i in x], squish_pts_means, width, label="SQUISH points kept (ratio 0.15)", color="#16a34a")
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel("Mean output points (avg over 10 large IDs)")
ax.set_title("Output size comparison (ratio 0.15 calibrated to SOTS)")
ax.legend()
ax.grid(True, axis="y", alpha=0.3)
for i, (s, q) in enumerate(zip(sots_pts_means, squish_pts_means)):
    ax.text(i - width/2, s+2, f"{s:.0f}", ha="center", va="bottom", fontsize=8, color="#2563eb")
    ax.text(i + width/2, q+2, f"{q:.0f}", ha="center", va="bottom", fontsize=8, color="#16a34a")

plt.tight_layout()
out4 = OUT_DIR / "points_bars.png"
plt.savefig(out4, dpi=180)
print(f"Wrote {out4}")
plt.close()

print("Done")
