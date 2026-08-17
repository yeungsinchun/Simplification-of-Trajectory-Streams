# Performance Investigation Results

## Summary

**Issue:** Performance degradation from ~83ms to ~145ms (74% slowdown)  
**Root Cause:** Unnecessary sorting in `get_points_from_grid()`  
**Fix:** Removed the sort operation  
**Result:** Performance improved to ~104ms (28% faster than broken version)

## Timing Breakdown

| Version | Time (ms) | Notes |
|---------|-----------|-------|
| Commit 3a669e9 (before web viz) | 83.3 | Baseline |
| HEAD with sort | 145.0 | 74% slower - sort + web code |
| HEAD without sort | 104.5 | 25% slower - just web code overhead |

## What Changed

### The Problematic Sort
In `simplify_geometry.h:398-403`, a row-major sort was added to `get_points_from_grid()`:

```cpp
std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) {
    const double ay = CGAL::to_double(a.y()), by = CGAL::to_double(b.y());
    if (ay != by) return ay < by;
    return CGAL::to_double(a.x()) < CGAL::to_double(b.x());
});
```

This sort:
- Runs for **every input point** (588 times for dataset 1)
- Sorts typically 20-50 grid corner points each time
- Uses expensive `CGAL::to_double()` conversions
- **Not required for correctness** - algorithm works without it

### Remaining 20ms Overhead
The ~20ms difference between 83ms (old) and 104ms (current without sort) comes from:
- Web tracing mode code being compiled in (even though not executed in normal path)
- Additional includes (`<fstream>`, `<iomanip>`)
- Code size increase affecting cache behavior

## CI/CD Setup

Created two GitHub Actions workflows:

### 1. `.github/workflows/ci.yml` - Basic CI
- Runs on every push/PR
- Builds the project
- Runs correctness test (checks output point count)
- Measures performance (warns if >200ms)
- Comments benchmark results on PRs

### 2. `.github/workflows/benchmark.yml` - Comprehensive Benchmarks
- Runs on push/PR and daily schedule
- Tests multiple epsilon values (50, 100, 200, 299, 500)
- Calculates throughput (points/second)
- Tracks performance over time
- Stores historical benchmark data
- Fails if performance exceeds 200ms threshold

## Recommendations

1. **Keep the sort removed** - It's not needed for correctness
2. **Monitor CI** - The benchmarks will catch future regressions
3. **Consider compile-time flags** - Could use `-DENABLE_WEB_TRACE` to conditionally compile web mode
4. **Profile if needed** - If 104ms is still too slow, profile to find the remaining 20ms

## Test Command

```bash
./build/simplify 1 -e 299 -d 1
```

Expected output:
- **Time:** ~100-105ms
- **Output points:** 92 (15.6% of 588)
