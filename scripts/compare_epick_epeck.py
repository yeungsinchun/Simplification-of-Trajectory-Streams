#!/usr/bin/env python3
import sys
import os

def load_points(path):
    pts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    pts.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    continue
    return pts

def frechet_distance(P, Q):
    """Standard DP Fréchet distance (free-space diagram)."""
    n, m = len(P), len(Q)
    if n == 0 or m == 0:
        return float('inf')

    # Pairwise squared distances
    d = [[0.0] * m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            dx = P[i][0] - Q[j][0]
            dy = P[i][1] - Q[j][1]
            d[i][j] = (dx*dx + dy*dy) ** 0.5

    ca = [[0.0] * m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            if i == 0 and j == 0:
                ca[i][j] = d[0][0]
            elif i > 0 and j == 0:
                ca[i][j] = max(ca[i-1][0], d[i][0])
            elif i == 0 and j > 0:
                ca[i][j] = max(ca[0][j-1], d[0][j])
            else:
                ca[i][j] = max(min(ca[i-1][j], ca[i][j-1], ca[i-1][j-1]), d[i][j])
    return ca[n-1][m-1]

data_dir = sys.argv[1] if len(sys.argv) > 1 else "/Users/sinchunyeung/Simplification-of-Trajectory-Streams/data/3"
epeck = load_points(os.path.join(data_dir, "simplify_epeck.txt"))
epick = load_points(os.path.join(data_dir, "simplify_epick.txt"))
original = load_points(os.path.join(data_dir, "original.txt"))

print(f"Original: {len(original)} points")
print(f"EPECK:    {len(epeck)} points")
print(f"EPICK:    {len(epick)} points")
print()

# Find differences
def find_diffs(a, b, tol=1e-8):
    diffs = []
    min_len = min(len(a), len(b))
    for i in range(min_len):
        dx = a[i][0] - b[i][0]
        dy = a[i][1] - b[i][1]
        dist = (dx*dx + dy*dy) ** 0.5
        if dist > tol:
            diffs.append((i, a[i], b[i], dist))
    return diffs, min_len

diffs, n = find_diffs(epeck, epick)
print(f"Points that differ (> 1e-8): {len(diffs)} out of {n}")
if diffs:
    print(f"\nFirst differing point at index {diffs[0][0]}:")
    print(f"  EPECK: ({diffs[0][1][0]:.12f}, {diffs[0][1][1]:.12f})")
    print(f"  EPICK: ({diffs[0][2][0]:.12f}, {diffs[0][2][1]:.12f})")
    print(f"  diff:  {diffs[0][3]:.12f}")
    print()
    for i, ep, ek, dist in diffs[:5]:
        print(f"  [{i:3d}] diff={dist:.10f}")
print()

# Check if they match at all
if len(epeck) == len(epick):
    all_match = all(
        ((e[0]-p[0])**2 + (e[1]-p[1])**2) ** 0.5 < 1e-8
        for e, p in zip(epeck, epick)
    )
    print(f"All {len(epeck)} points match (tol=1e-8): {all_match}")

    if all_match:
        print("\nEPICK and EPECK produce IDENTICAL outputs on this dataset.")
        print("The difference in reported Fréchet must come from how simplify.cpp")
        print("computes the distance (different code path, not just the kernel).")
        print()

# Fréchet distance to original
print("Fréchet distance to original trajectory:")
fd_orig_epeck = frechet_distance(original, epeck)
fd_orig_epick = frechet_distance(original, epick)
print(f"  Fréchet(Original, EPECK_simplified) = {fd_orig_epeck:.6f}")
print(f"  Fréchet(Original, EPICK_simplified) = {fd_orig_epick:.6f}")
print(f"  Difference = {abs(fd_orig_epeck - fd_orig_epick):.6f}")

# Check the first and last few simplified points
print("\nFirst 3 simplified points:")
for i in range(min(3, len(epeck), len(epick))):
    print(f"  [{i}] EPECK=({epeck[i][0]:.4f}, {epeck[i][1]:.4f})  EPICK=({epick[i][0]:.4f}, {epick[i][1]:.4f})")

print("\nLast 3 simplified points:")
for offset in range(1, 4):
    i = len(epeck) - offset
    j = len(epick) - offset
    if i >= 0 and j >= 0:
        print(f"  [{i}] EPECK=({epeck[i][0]:.4f}, {epeck[i][1]:.4f})  EPICK=({epick[j][0]:.4f}, {epick[j][1]:.4f})")
