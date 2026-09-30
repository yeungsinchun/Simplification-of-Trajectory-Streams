#!/usr/bin/env bash
# Decide whether a change could plausibly affect benchmark performance.
#
# Usage: perf_paths.sh <base-sha> <head-sha>
# Prints `perf=true|false` (suitable for appending to $GITHUB_OUTPUT) and
# the matching paths on stderr. Fails closed: an unusable diff yields true.
#
# Keep this list explicit. A path belongs here only if it feeds the `simplify`
# binary, its build, the benchmark inputs, or the benchmark/gating logic.
set -uo pipefail

is_perf_path() {
  case "$1" in
    # Core algorithm sources linked into `simplify` (simplify.cpp,
    # web_trace.cpp, simplify_*.h, timer.h, web_trace.h).
    *.cpp | *.h) [[ "$1" != */* ]] ;;
    # Build definition, toolchain flags and the baseline-binary submodule.
    CMakeLists.txt | traj-compression | traj-compression/*) return 0 ;;
    # Benchmark trajectories and their derivation.
    data/* | scripts/derive_benchmark_data.py) return 0 ;;
    # Benchmark statistics, reporting and their tests.
    scripts/ci/welch.py | scripts/ci/bench_report.py | scripts/ci/test_bench_report.py) return 0 ;;
    # The benchmark workflow itself and this gate.
    .github/workflows/benchmark.yml | scripts/ci/perf_paths.sh) return 0 ;;
    *) return 1 ;;
  esac
}

base="${1:-}"
head="${2:-}"

emit() { echo "perf=$1"; }

if [ -z "$base" ] || [ -z "$head" ] || [[ "$base" =~ ^0+$ ]]; then
  echo "no usable base commit; running benchmark" >&2
  emit true
  exit 0
fi

if ! files=$(git diff --name-only "$base" "$head" --); then
  echo "git diff failed; running benchmark" >&2
  emit true
  exit 0
fi

perf=false
while IFS= read -r f; do
  [ -n "$f" ] || continue
  if is_perf_path "$f"; then
    echo "performance-relevant: $f" >&2
    perf=true
  fi
done <<< "$files"

[ "$perf" = true ] || echo "no performance-relevant paths changed" >&2
emit "$perf"
