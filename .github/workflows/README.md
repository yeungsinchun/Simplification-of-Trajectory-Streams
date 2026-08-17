# GitHub Actions Workflows

This directory contains CI/CD workflows for the trajectory simplification project.

## Workflows

### ci.yml - Basic CI
**Triggers:** Push to any branch, Pull requests
**Duration:** ~2-3 minutes

Tests:
- Builds the project in Release mode
- Runs correctness test on dataset 1
- Measures performance baseline (ε=299, δ=1)
- Posts results as PR comment

Pass criteria:
- Build succeeds
- Output has correct number of points (92 for dataset 1)
- Performance is under 200ms threshold

### benchmark.yml - Comprehensive Benchmarks
**Triggers:** Push to main, Pull requests, Daily at 00:00 UTC, Manual dispatch
**Duration:** ~3-5 minutes

Tests:
- Builds in Release mode
- Tests 5 different epsilon values (50, 100, 200, 299, 500)
- Calculates compression ratios and throughput
- Stores historical data on main branch
- Posts detailed report as PR comment

Pass criteria:
- All builds succeed
- Baseline performance (ε=299) under 200ms
- Results stored in benchmark_history/ for trending

## Expected Performance

With Release build on modern hardware:
- **Target:** ~10-20ms for 588 points
- **Threshold:** 200ms (workflow fails if exceeded)
- **Throughput:** ~30,000-60,000 points/second

## Local Testing

To run the same tests locally:

```bash
# Build
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --target simplify -j$(nproc)

# Test
cd ..
./build/simplify 1 -e 299 -d 1

# Should see:
# Loaded 588 points.
# SIMPLIFY_CORE_MS: 10-20
# Simplified to 92 points (15.6%)
```
