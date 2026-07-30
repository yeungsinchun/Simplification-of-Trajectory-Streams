# Algorithmic Complexity Analysis: Simplify vs DOTS

## Executive Summary

**Your implementation is ~700x slower than DOTS** because it uses fundamentally different algorithms with vastly different computational complexity.

---

## Operation Count Comparison

### DOTS Algorithm (Fast - LSSD Metric)

**Per-point operations:**
- **Line-Segment-to-Segment Distance (LSSD)**: O(1) geometric distance calculation
- **Simple comparison**: Check if LSSD exceeds threshold

**Total complexity for trajectory of length n:**
```
O(n) time complexity
O(1) space complexity per iteration
```

**Concrete operations for n=10,000 points:**
- ~10,000 LSSD distance calculations
- ~10,000 comparisons
- **Estimated: ~20,000 simple operations**

---

### Your Implementation (Slow - Geometric Intersection)

**Per-point operations (for each candidate in grid):**

1. **Grid generation** - `get_points_from_grid()`:
   - Grid size: ~(2R/GRID)² cells
   - With δ=200, ε=0.5: R ≈ 250, GRID ≈ 35
   - Result: ~50-100 grid corners per point
   - Complexity: **O(R²/GRID²)** per point

2. **Convex hull** - `get_conv_from_grid()`:
   - CGAL convex hull on ~50-100 points
   - Complexity: **O(k log k)** where k ≈ 50-100

3. **For EACH of ~50-100 grid candidates, EACH input point:**
   
   a. **find_F()** - Free-space wedge construction:
      - Find tangent vertices: O(n) scan through stab region
      - Ray-bbox intersection: O(1)
      - Wedge stitching: O(n)
      - **Total: O(n)** per candidate
   
   b. **intersect()** - Sutherland-Hodgman polygon clipping:
      - `dedup_consecutive()`: O(n + m)
      - `clip()`: **O(n × m)** where n,m = polygon vertex counts
        - Inner loop: for each of m clip edges, process n subject vertices
        - With n,m ≈ 10-20: ~200-400 operations per clip
      - **Total: O(n × m)** per intersection

**Total per input point consumed:**
```
For Pn candidates (typically 50-100):
  Pn × [find_F: O(S) + intersect: O(F × G)]
  
Where:
  S = stab region size (grows with trajectory, can be 10-50 vertices)
  F = wedge polygon size (~10-20 vertices)
  G = grid convex hull size (~8-16 vertices)
```

**Conservative estimate for one input point:**
- 50 candidates × (20 ops for find_F + 200 ops for clip) = **11,000 operations**

**For n=10,000 points:**
- 10,000 points × 11,000 ops/point = **110,000,000 operations**
- Compare to DOTS: 20,000 operations
- **Ratio: 110M / 20K = 5,500x** (theoretical worst case)

---

## Why the 700x Slowdown?

### 1. **Algorithmic Complexity Gap**

| Metric | DOTS | Your Implementation |
|--------|------|---------------------|
| Per-point cost | O(1) | O(Pn × S × F × G) |
| Typical operations | 2 | ~11,000 |
| Memory allocations | 0 | 4-6 per point |

### 2. **Polygon Clipping Overhead**

The Sutherland-Hodgman `clip()` function:
- Called **~320,000 times** for test case 1 (n≈6,400 points × 50 candidates)
- Each call: 4 vector allocations + O(n×m) loop iterations
- **Before optimization**: Each call allocates `subj`, `clipv`, `out`, `in`, `result`
- **After optimization**: Reuses workspace, saves ~1.6M allocations

### 3. **Geometric Precision vs Speed Trade-off**

**DOTS**: 
- Simple Euclidean distance in doubles
- No geometric primitives, no CGAL overhead

**Your implementation**:
- Exact predicates (CGAL kernel)
- Convex hull computations
- Polygon intersection with deduplication
- Ray-tracing, tangent finding, wedge construction

---

## Benchmark Results Explained

### Test Case 1 (n=6,400):
```
DOTS:     9.10 ms
Simplify: 6,409.52 ms
Slowdown: 704x
```

**Why not 5,500x?**
- Many candidates die early (dead_cnt optimization)
- Polygons are small in practice (F,G ≈ 8-12 vertices, not 20)
- CGAL is highly optimized
- **Still fundamentally O(n × k × poly_ops) vs O(n)**

### Test Case 100 (different trajectory):
```
DOTS:     8.82 ms
Simplify: 6,043.46 ms
Slowdown: 685x
```

Consistent ~700x ratio confirms this is **algorithmic**, not implementation fluff.

---

## Performance Improvements

### ✅ Already Implemented: Workspace Reuse

**Before:**
```cpp
inline std::vector<Point> clip(...) {
    std::vector<std::array<double, 2>> subj, clipv;  // allocate
    std::vector<std::array<double, 2>> out = subj;    // allocate
    std::vector<std::array<double, 2>> in;            // allocate
    std::vector<Point> result;                        // allocate
    // ... computation ...
}
```

**After:**
```cpp
thread_local ClipWorkspace ws;  // ONE allocation per thread, reused forever

inline std::vector<Point> clip(...) {
    ClipWorkspace& ws = get_clip_workspace();
    ws.subj.clear(); ws.clipv.clear(); // reuse existing capacity
    // ... computation ...
}
```

**Expected speedup:** 30-50% reduction in time (eliminates ~1.6M allocations)

---

### 🎯 Other Potential Optimizations

1. **Early bounding-box rejection** (5-10% speedup):
   ```cpp
   if (!bbox_overlap(F, Gi)) continue;  // skip expensive clip()
   ```

2. **Batch-kill dead candidates** (marginal):
   - Current: check `if (dead[i])` every iteration
   - Better: maintain active_list, remove dead ones

3. **Spatial indexing** (complex, 2-3x speedup):
   - Use R-tree or grid to avoid checking all Pn candidates

4. **SIMD vectorization** (expert-level, 2-4x speedup):
   - Vectorize the inner loops of Sutherland-Hodgman

---

## Fundamental Reality

**You cannot close the 700x gap without changing the algorithm.**

The geometric intersection approach is solving a **harder problem** than DOTS:
- DOTS: "Is this point too far from the line segment?"
- Yours: "What is the exact geometric region reachable through this wedge?"

Your algorithm provides **stronger guarantees** (formal Fréchet distance bound), but at **massive computational cost**.

---

## Recommendations

1. ✅ **Use the workspace optimization** (30-50% speedup expected)
2. ✅ **Profile with xcrun xctrace** to find remaining hot spots
3. ⚠️ **Accept the performance gap** if geometric correctness is required
4. ❌ **Don't chase micro-optimizations** - you're fighting O(n²) vs O(n)

If speed matters more than geometric precision, **use DOTS or implement LSSD-based simplification**.
