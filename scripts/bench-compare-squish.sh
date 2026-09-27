#!/usr/bin/env bash
# bench-compare-squish.sh — Streaming SOTS vs SQUISH baseline benchmark.
# Builds both binaries, runs 10-sample Welford means per epsilon tier,
# computes speedup vs SQUISH with Welch t, writes JSON/CSV to .lavish.
#
# Algorithmic complexity (papers/journal.pdf, papers/squish.pdf):
#   SOTS: O(epsilon^{-alpha}) storage, O(epsilon^{-alpha} log 1/epsilon) per vertex for d=2,3
#         where alpha=2(d-1) floor(d/2)^2+d (d=2 => alpha=4 => O(epsilon^{-4})).
#         Guarantee dF <= (1+epsilon)delta and |sigma| <=2 kappa(delta)-2. Poly(1/epsilon) geometry.
#   SQUISH: streaming heuristic, buffer B=ratio*N, O(B) per point naive (O(log B) heap), O(N*B) total, O(B) storage, SED, no Frechet bound.
#   DP: offline batch, O(N log N) avg O(N^2) worst, PED, no guarantee. DOTS: O(N/M) per point.
#   SOTS pays poly(1/epsilon) for guarantee; extra-fine epsilon=0.1 is 10,000x geometry of epsilon=299.
#
# Parameter tuning (fair before results):
#   SOTS: epsilon in {299,30,5,0.5,0.1} (5 tiers extra-coarse..extra-fine), delta=NUM/(1+epsilon) NUM=300 constant,
#         so (1+epsilon)delta=300 fixed — isolates approximation tightness from envelope size.
#         Smaller epsilon => finer grid epsilon*delta/(2 sqrt(d)) => tighter (1+epsilon)delta bound but slower.
#   SQUISH: ratio sweep 0.05/0.15/0.5 on data/21 (N=588). SOTS extra-coarse keeps ~15% (86-92 pts).
#           ratio 0.15 keeps ~15% (88 pts) — size-matched fair. 0.5 keeps 3x more (not fair). 0.15 is default.
#   Measurement: same IDs, same order, CORE_MS (no I/O), 10-run Welford, Welch t 95% — same as benchmark.yml.
#
# Usage: scripts/bench-compare-squish.sh [--epsilon "0.1 0.5 ..."] [--delta-numer 300] [--size large] [--runs 10] [--ratio 0.15]
# Any contributor can run: ./scripts/bench-compare-squish.sh
# Output: .lavish/bench-squish-*.json, .lavish/bench-squish-*.csv and markdown table on stdout.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${REPO_ROOT}/build"
LAVISH_DIR="${REPO_ROOT}/.lavish"
WELCH_PY="${REPO_ROOT}/scripts/ci/welch.py"

# Defaults mirror CI matrix: 5 epsilon tiers, delta 300/(1+e), large size.
EPSILONS="299 30 5 0.5 0.1"
DELTA_NUMER="300"
SIZE="large"
IDS_SMALL="11 12 13 14 15 16 17 18 19 20"
IDS_LARGE="21 22 23 24 25 26 27 28 29 30"
RUNS=10
RATIO="0.15"
OUTPUT_DIR=""
WARMUP=1

print_help() {
  cat <<'HELP'
Usage: scripts/bench-compare-squish.sh [options]

Benchmark streaming SOTS (simplify) vs SQUISH baseline at multiple
epsilon tiers. Reuses CI gating methodology: 10-run Welford mean/stddev,
Welch t confidence, same-binary comparison disciplines.

Complexity (papers/journal.pdf):
  SOTS  O(epsilon^{-alpha}) storage, O(epsilon^{-alpha} log 1/epsilon) per vertex (d=2,3),
        alpha=2(d-1) floor(d/2)^2+d, d=2 => alpha=4 => O(epsilon^{-4}), guarantee dF<=(1+epsilon)delta.
  SQUISH O(B) per point naive B=ratio*N, O(N*B) total, O(B) storage, SED heuristic, no bound.
  DP     O(N log N) avg O(N^2) worst, batch, no bound.  SOTS trades poly(1/epsilon) for guarantee.

Parameter tuning (before results):
  SOTS epsilon {299,30,5,0.5,0.1}, delta=NUM/(1+epsilon) NUM=300 so (1+epsilon)delta constant=300.
       Smaller epsilon => finer grid epsilon*delta/(2 sqrt(d)) => tighter bound but slower (extra-fine 10,000x geometry).
  SQUISH ratio 0.15 calibrated: swept 0.05/0.15/0.5 on data/21 (N=588), SOTS extra-coarse ~15% (86-92 pts),
       SQUISH 0.15 => 88 pts (size-matched, fair), 0.5 => 294 pts (3x more, not fair). Fixed before timing.

Options:
  --epsilon "LIST"     Space-separated epsilon tiers (default: "299 30 5 0.5 0.1")
  --delta-numer NUM    Numerator for delta = NUM/(1+epsilon) (default: 300)
  --delta-numer-list "LIST"  Space-separated numerators (default: 300, use "300 1000" for full matrix)
  --size SIZE          small, large, or all (default: large)
  --ids "LIST"         Override ID list (default: derived from --size)
  --runs N             Samples per ID/binary (default: 10)
  --ratio R            SQUISH keep ratio 0..1 (default: 0.15, matches SOTS output at extra-coarse)
  --output-dir DIR     Override output directory (default: .lavish)
  --no-warmup          Skip warmup run
  -h, --help           Show this help

Examples:
  ./scripts/bench-compare-squish.sh
  ./scripts/bench-compare-squish.sh --epsilon "0.1 0.5" --ratio 0.5
  ./scripts/bench-compare-squish.sh --epsilon "299" --runs 5 --size small

Output:
  .lavish/bench-squish-YYYYmmdd-HHMMSS.json
  .lavish/bench-squish-YYYYmmdd-HHMMSS.csv
  Markdown table printed to stdout.

Repro: uses SIMPLIFY_CORE_MS / SQUISH_CORE_MS (core algorithm wall ms, no I/O)
       exactly as benchmark.yml does, with Welford + Welch via scripts/ci/welch.py.

HELP
}

DELTA_NUMERS=""

# Parse args
while [ $# -gt 0 ]; do
  case "$1" in
    --epsilon)
      EPSILONS="${2:-}"
      shift 2
      ;;
    --delta-numer)
      DELTA_NUMER="${2:-}"
      DELTA_NUMERS="${2:-}"
      shift 2
      ;;
    --delta-numer-list)
      DELTA_NUMERS="${2:-}"
      shift 2
      ;;
    --size)
      SIZE="${2:-}"
      shift 2
      ;;
    --ids)
      IDS_OVERRIDE="${2:-}"
      shift 2
      ;;
    --runs)
      RUNS="${2:-}"
      shift 2
      ;;
    --ratio)
      RATIO="${2:-}"
      shift 2
      ;;
    --output-dir)
      OUTPUT_DIR="${2:-}"
      shift 2
      ;;
    --no-warmup)
      WARMUP=0
      shift
      ;;
    -h|--help)
      print_help
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      print_help >&2
      exit 2
      ;;
  esac
done

# Validate
if ! echo "$RUNS" | grep -Eq '^[0-9]+$' || [ "$RUNS" -lt 2 ]; then
  echo "error: --runs must be integer >=2" >&2; exit 2
fi
if ! python3 -c "r=float('${RATIO}'); assert 0 < r <= 1" 2>/dev/null; then
  echo "error: --ratio must be 0 < ratio <= 1" >&2; exit 2
fi

# Determine IDs based on size if not overridden
if [ -n "${IDS_OVERRIDE:-}" ]; then
  IDS="$IDS_OVERRIDE"
else
  case "$SIZE" in
    small) IDS="$IDS_SMALL" ;;
    large) IDS="$IDS_LARGE" ;;
    all)   IDS="$IDS_SMALL $IDS_LARGE" ;;
    *) echo "error: --size must be small|large|all" >&2; exit 2 ;;
  esac
fi

# Delta numerators handling
if [ -z "$DELTA_NUMERS" ]; then
  DELTA_NUMERS="$DELTA_NUMER"
fi

# Output dir
if [ -n "$OUTPUT_DIR" ]; then
  LAVISH_DIR="$OUTPUT_DIR"
fi
mkdir -p "$LAVISH_DIR"

# Helpers
compute_delta() {
  # $1=epsilon $2=numer
  python3 -c "import sys; print('{:.15g}'.format(float(sys.argv[2])/(1+float(sys.argv[1]))))" "$1" "$2"
}

bench_stats_sots() {
  # $1=id $2=epsilon $3=delta $4=runs -> prints "mean stddev"
  local id="$1" eps="$2" delta="$3" n="$4"
  local count=0 mean=0 m2=0 ms run log
  for run in $(seq 1 "$n"); do
    log="$(mktemp)"
    if ! ( cd "${REPO_ROOT}" && ./build/simplify "${id}" -e "${eps}" -d "${delta}" >"${log}" 2>&1 ); then
      echo "❌ SOTS run failed id=${id} e=${eps} d=${delta} run ${run}/${n}" >&2
      cat "${log}" >&2
      rm -f "${log}"
      return 1
    fi
    ms="$(grep -m1 '^SIMPLIFY_CORE_MS:' "${log}" | awk '{print $2}')"
    rm -f "${log}"
    if [ -z "${ms}" ]; then
      echo "❌ Missing SIMPLIFY_CORE_MS id=${id} run ${run}/${n}" >&2
      return 1
    fi
    count=$((count + 1))
    # Welford
    delta_w=$(awk -v x="${ms}" -v m="${mean}" 'BEGIN{printf "%.17g", x - m}')
    mean=$(awk -v m="${mean}" -v d="${delta_w}" -v c="${count}" 'BEGIN{printf "%.17g", m + d / c}')
    delta_w2=$(awk -v x="${ms}" -v m="${mean}" 'BEGIN{printf "%.17g", x - m}')
    m2=$(awk -v cur="${m2}" -v d="${delta_w}" -v d2="${delta_w2}" 'BEGIN{printf "%.17g", cur + d * d2}')
  done
  local variance=0 stddev=0
  if [ "${count}" -gt 1 ]; then
    variance=$(awk -v m2="${m2}" -v c="${count}" 'BEGIN{printf "%.17g", m2/(c-1)}')
    variance=$(awk -v v="${variance}" 'BEGIN{if(v<0)v=0; printf "%.17g", v}')
    stddev=$(awk -v v="${variance}" 'BEGIN{printf "%.17g", sqrt(v)}')
  fi
  printf '%s %s' "${mean}" "${stddev}"
}

bench_stats_squish() {
  # $1=id $2=ratio $3=runs -> prints "mean stddev"
  local id="$1" ratio="$2" n="$3"
  local count=0 mean=0 m2=0 ms run log
  for run in $(seq 1 "$n"); do
    log="$(mktemp)"
    out="$(mktemp)"
    if ! ( ./build/squish "data/${id}/original.txt" "${ratio}" "${out}" >"${log}" 2>&1 ); then
      echo "❌ SQUISH run failed id=${id} ratio=${ratio} run ${run}/${n}" >&2
      cat "${log}" >&2
      rm -f "${log}" "${out}"
      return 1
    fi
    ms="$(grep -m1 '^SQUISH_CORE_MS:' "${log}" | awk '{print $2}')"
    rm -f "${log}" "${out}"
    if [ -z "${ms}" ]; then
      # Fallback: try wall time if CORE missing (should not happen)
      echo "❌ Missing SQUISH_CORE_MS id=${id} run ${run}/${n}" >&2
      return 1
    fi
    count=$((count + 1))
    delta_w=$(awk -v x="${ms}" -v m="${mean}" 'BEGIN{printf "%.17g", x - m}')
    mean=$(awk -v m="${mean}" -v d="${delta_w}" -v c="${count}" 'BEGIN{printf "%.17g", m + d / c}')
    delta_w2=$(awk -v x="${ms}" -v m="${mean}" 'BEGIN{printf "%.17g", x - m}')
    m2=$(awk -v cur="${m2}" -v d="${delta_w}" -v d2="${delta_w2}" 'BEGIN{printf "%.17g", cur + d * d2}')
  done
  local variance=0 stddev=0
  if [ "${count}" -gt 1 ]; then
    variance=$(awk -v m2="${m2}" -v c="${count}" 'BEGIN{printf "%.17g", m2/(c-1)}')
    variance=$(awk -v v="${variance}" 'BEGIN{if(v<0)v=0; printf "%.17g", v}')
    stddev=$(awk -v v="${variance}" 'BEGIN{printf "%.17g", sqrt(v)}')
  fi
  printf '%s %s' "${mean}" "${stddev}"
}

# Build check
echo "Building simplify and squish (Release, BUILD_GUI=OFF)..." >&2
if [ ! -f "${REPO_ROOT}/build/simplify" ] || [ ! -f "${REPO_ROOT}/build/squish" ]; then
  cmake -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release -DBUILD_GUI=OFF > /dev/null
  cmake --build "${BUILD_DIR}" --target simplify squish -j "$( (nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4) )" > /dev/null
fi
# Rebuild to ensure latest
cmake -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release -DBUILD_GUI=OFF > /dev/null 2>&1 || true
cmake --build "${BUILD_DIR}" --target simplify squish -j "$( (nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4) )" > /dev/null 2>&1 || {
  echo "Build failed, trying verbose..." >&2
  cmake --build "${BUILD_DIR}" --target simplify squish -j 4
}

# Verify data exists
missing_ids=""
for id in ${IDS}; do
  if [ ! -f "${REPO_ROOT}/data/${id}/original.txt" ]; then
    missing_ids="${missing_ids} ${id}"
  fi
done
if [ -n "${missing_ids}" ]; then
  echo "error: missing data for IDs:${missing_ids}. Run python3 scripts/derive_benchmark_data.py" >&2
  exit 1
fi

# Hardware info for JSON
HOST_INFO="$(uname -a 2>/dev/null || echo unknown)"
COMPILER_INFO="$({ c++ --version 2>&1 | head -n1; } || echo unknown)"
CGAL_INFO="$({ grep -r "CGAL.*version" /opt/homebrew/include/CGAL/version.h 2>/dev/null | head -n1; } || echo "CGAL latest")"
if [ -f /opt/homebrew/include/CGAL/version.h ]; then
  CGAL_INFO="$(grep -E "#define CGAL_VERSION" /opt/homebrew/include/CGAL/version.h | head -n1 || echo CGAL)"
fi
GIT_SHA="$(git -C "${REPO_ROOT}" rev-parse --short HEAD 2>/dev/null || echo unknown)"
DATE_ISO="$(date -u +%Y-%m-%dT%H:%M:%SZ 2>/dev/null || date +%Y%m%d-%H%M%S)"
STAMP="$(date +%Y%m%d-%H%M%S 2>/dev/null || date -u +%Y%m%d%H%M%S)"
JSON_OUT="${LAVISH_DIR}/bench-squish-${STAMP}.json"
CSV_OUT="${LAVISH_DIR}/bench-squish-${STAMP}.csv"
TSV_OUT="${LAVISH_DIR}/bench-squish-${STAMP}.tsv"

# Warmup
if [ "${WARMUP}" -eq 1 ]; then
  echo "Warming up (1 run per binary)..." >&2
  first_id="$(echo "${IDS}" | awk '{print $1}')"
  first_eps="$(echo "${EPSILONS}" | awk '{print $1}')"
  first_delta="$(compute_delta "${first_eps}" "$(echo "${DELTA_NUMERS}" | awk '{print $1}')")"
  ./build/simplify "${first_id}" -e "${first_eps}" -d "${first_delta}" > /dev/null 2>&1 || true
  ./build/squish "data/${first_id}/original.txt" "${RATIO}" "/tmp/warm_out.txt" > /dev/null 2>&1 || true
  rm -f /tmp/warm_out.txt
fi

# Prepare CSV header
printf 'epsilon,delta,delta_numer,id,sots_ms,sots_std,squish_ms,squish_std,ratio,speedup_squish_over_sots,sots_points,squish_points\n' > "${CSV_OUT}"
printf 'epsilon\tdelta\tdelta_numer\tid\tsots_ms\tsots_std\tsquish_ms\tsquish_std\tratio\tspeedup\tsots_points\tsquish_points\n' > "${TSV_OUT}"

# Collect for markdown
MARKDOWN_ROWS=""
JSON_CASES_TMP="$(mktemp)"
echo "[]" > "${JSON_CASES_TMP}"

# overall aggregates collected via markdown rows and JSON

# Loop over delta numerators and epsilons
for delta_numer in ${DELTA_NUMERS}; do
  for eps in ${EPSILONS}; do
    delta="$(compute_delta "${eps}" "${delta_numer}")"
    label=""
    # friendly label like extra-fine-e-d300-large
    case "${eps}" in
      299) tier="extra-coarse-e" ;;
      30) tier="coarse-e" ;;
      5) tier="mid-e" ;;
      0.5) tier="fine-e" ;;
      0.1) tier="extra-fine-e" ;;
      *) tier="e${eps}" ;;
    esac
    label="${tier}-d${delta_numer}-${SIZE}"

    echo "Benchmarking ${label}: epsilon=${eps} delta=${delta} (${delta_numer}/(1+${eps})) runs=${RUNS} ids=${IDS} ratio=${RATIO}" >&2

    per_eps_sots_means=()
    per_eps_sots_stds=()
    per_eps_squish_means=()
    per_eps_squish_stds=()
    per_eps_sots_pts=()
    per_eps_squish_pts=()

    for id in ${IDS}; do
      # SOTS
      if ! sots_stats="$(bench_stats_sots "${id}" "${eps}" "${delta}" "${RUNS}")"; then
        echo "Failed SOTS bench id=${id} eps=${eps}" >&2
        exit 1
      fi
      sots_ms="$(echo "${sots_stats}" | awk '{print $1}')"
      sots_std="$(echo "${sots_stats}" | awk '{print $2}')"
      # SQUISH
      if ! squish_stats="$(bench_stats_squish "${id}" "${RATIO}" "${RUNS}")"; then
        echo "Failed SQUISH bench id=${id}" >&2
        exit 1
      fi
      squish_ms="$(echo "${squish_stats}" | awk '{print $1}')"
      squish_std="$(echo "${squish_stats}" | awk '{print $2}')"

      # Points: run once to get output sizes
      sots_out="/tmp/sots_${id}_pts.txt"
      squish_out="/tmp/squish_${id}_pts.txt"
      # Use simplify to produce output (we can just count after run)
      # Simplify writes to data/<id>/simplify.txt, read header
      # Run simplify once more to ensure output corresponds to this epsilon
      ./build/simplify "${id}" -e "${eps}" -d "${delta}" > /dev/null 2>&1 || true
      sots_pts="0"
      if [ -f "data/${id}/simplify.txt" ]; then
        sots_pts="$(head -n1 "data/${id}/simplify.txt" 2>/dev/null | tr -d '\r' || echo 0)"
      fi
      ./build/squish "data/${id}/original.txt" "${RATIO}" "${squish_out}" > /dev/null 2>&1 || true
      squish_pts="$(head -n1 "${squish_out}" 2>/dev/null | tr -d '\r' || echo 0)"
      rm -f "${sots_out}" "${squish_out}"

      # Ratio and speedup: squish / sots (how many times faster squish is)
      # Use python for precise
      ratio_val="$(python3 -c "import sys; s=float(sys.argv[1]); q=float(sys.argv[2]); print('{:.6g}'.format(q/s if s!=0 else float('inf')))" "${sots_ms}" "${squish_ms}")"
      speedup_val="${ratio_val}"
      # For markdown, show sots vs squish
      # Also compute inverse for table if needed

      printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' "${eps}" "${delta}" "${delta_numer}" "${id}" "${sots_ms}" "${sots_std}" "${squish_ms}" "${squish_std}" "${RATIO}" "${speedup_val}" "${sots_pts}" "${squish_pts}" >> "${CSV_OUT}"
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${eps}" "${delta}" "${delta_numer}" "${id}" "${sots_ms}" "${sots_std}" "${squish_ms}" "${squish_std}" "${RATIO}" "${speedup_val}" "${sots_pts}" "${squish_pts}" >> "${TSV_OUT}"

      per_eps_sots_means+=("${sots_ms}")
      per_eps_sots_stds+=("${sots_std}")
      per_eps_squish_means+=("${squish_ms}")
      per_eps_squish_stds+=("${squish_std}")
      per_eps_sots_pts+=("${sots_pts}")
      per_eps_squish_pts+=("${squish_pts}")

      # JSON per case via python
      python3 - "${JSON_CASES_TMP}" "${eps}" "${delta}" "${delta_numer}" "${id}" "${sots_ms}" "${sots_std}" "${squish_ms}" "${squish_std}" "${RATIO}" "${speedup_val}" "${sots_pts}" "${squish_pts}" <<'PY'
import json, sys
tmp, eps, delta, numer, id_, sots_ms, sots_std, squish_ms, squish_std, ratio, speedup, sots_pts, squish_pts = sys.argv[1:]
with open(tmp) as f:
    arr = json.load(f)
arr.append({
    "label": f"e{eps}-d{numer}",
    "epsilon": float(eps),
    "delta": float(delta),
    "delta_numer": float(numer),
    "id": int(id_),
    "sots_ms": float(sots_ms),
    "sots_std": float(sots_std),
    "squish_ms": float(squish_ms),
    "squish_std": float(squish_std),
    "ratio": float(ratio),
    "speedup_squish_over_sots": float(speedup) if speedup != "inf" else None,
    "sots_points": int(sots_pts) if sots_pts.isdigit() else None,
    "squish_points": int(squish_pts) if squish_pts.isdigit() else None,
})
with open(tmp, "w") as f:
    json.dump(arr, f, indent=2)
PY

    done

    # Compute per-epsilon aggregate means
    # Mean sots, squish across ids
    mean_sots="$(python3 -c "import sys; vals=list(map(float, sys.argv[1:])); print(sum(vals)/len(vals) if vals else 0)" "${per_eps_sots_means[@]}")"
    mean_squish="$(python3 -c "import sys; vals=list(map(float, sys.argv[1:])); print(sum(vals)/len(vals) if vals else 0)" "${per_eps_squish_means[@]}")"
    mean_speedup="$(python3 -c "import sys; a=float(sys.argv[1]); b=float(sys.argv[2]); print(b/a if a!=0 else 0)" "${mean_sots}" "${mean_squish}")"
    mean_sots_str="$(printf '%.4f' "${mean_sots}")"
    mean_squish_str="$(printf '%.4f' "${mean_squish}")"
    mean_speedup_str="$(python3 -c "import sys; print('{:.3f}'.format(float(sys.argv[1])))" "${mean_speedup}")"

    # Welch verdict: is squish significantly faster? Use welch.py with limit 1.0
    # verdict = 1 if (squish - limit*sots) lower bound >0 => squish slower, 0 otherwise
    # We want p-like: test sots vs squish. For now compute welch for sots vs squish with limit 1.0 both ways.
    # Use mean comparison: compare squish vs sots*1.0
    welch_speedup_p="n/a"
    welch_verdict="n/a"
    if [ -f "${WELCH_PY}" ]; then
      # Prepare space-separated lists
      sots_means_str="${per_eps_sots_means[*]}"
      sots_stds_str="${per_eps_sots_stds[*]}"
      squish_means_str="${per_eps_squish_means[*]}"
      squish_stds_str="${per_eps_squish_stds[*]}"
      # Welch: is squish > sots? (i.e., squish slower) test with limit 1.0
      # We run: welch.py RUNS 1.0 "sots_means" "sots_stds" "squish_means" "squish_stds"
      # If verdict 1, squish is confidently slower (mean_squish > mean_sots)
      if welch_out="$(python3 "${WELCH_PY}" "${RUNS}" "1.0" "${sots_means_str}" "${sots_stds_str}" "${squish_means_str}" "${squish_stds_str}" 2>&1)"; then
        welch_verdict="${welch_out}"
        if [ "${welch_verdict}" = "1" ]; then
          welch_speedup_p="squish slower (p<0.05)"
        elif [ "${welch_verdict}" = "0" ]; then
          # Also test opposite: sots > squish ?
          if welch_out2="$(python3 "${WELCH_PY}" "${RUNS}" "1.0" "${squish_means_str}" "${squish_stds_str}" "${sots_means_str}" "${sots_stds_str}" 2>&1)"; then
            if [ "${welch_out2}" = "1" ]; then
              welch_speedup_p="sots slower (p<0.05)"
            else
              welch_speedup_p="no significant diff"
            fi
          fi
        fi
      else
        welch_speedup_p="error"
      fi
    fi

    # For markdown rows, compute per-epsilon std aggregated? Use mean std as avg std
    avg_sots_std="$(python3 -c "import sys; vals=list(map(float, sys.argv[1:])); print(sum(vals)/len(vals) if vals else 0)" "${per_eps_sots_stds[@]}")"
    avg_squish_std="$(python3 -c "import sys; vals=list(map(float, sys.argv[1:])); print(sum(vals)/len(vals) if vals else 0)" "${per_eps_squish_stds[@]}")"
    printf -v avg_sots_std_str '%.4f' "${avg_sots_std}"
    printf -v avg_squish_std_str '%.4f' "${avg_squish_std}"

    # Store for markdown
    MARKDOWN_ROWS+="${eps}|${delta}|${delta_numer}|${mean_sots_str} ± ${avg_sots_std_str}|${mean_squish_str} ± ${avg_squish_std_str}|${mean_speedup_str}×|${welch_speedup_p}
"
    # Also accumulate for JSON config
    python3 - "${JSON_CASES_TMP}" "${eps}" "${delta}" "${delta_numer}" "${label}" "${mean_sots}" "${mean_squish}" "${mean_speedup}" "${welch_verdict}" <<'PY'
import json, sys
tmp, eps, delta, numer, label, mean_sots, mean_squish, mean_speedup, welch = sys.argv[1:]
with open(tmp) as f:
    data = json.load(f)
# Append aggregate as separate entry with id=-1
data.append({
    "_aggregate": True,
    "label": label,
    "epsilon": float(eps),
    "delta": float(delta),
    "delta_numer": float(numer),
    "mean_sots_ms": float(mean_sots),
    "mean_squish_ms": float(mean_squish),
    "mean_speedup_squish_over_sots": float(mean_speedup),
    "welch_verdict": welch,
})
with open(tmp, "w") as f:
    json.dump(data, f, indent=2)
PY

  done
done

# Build final JSON with meta
python3 - "${JSON_CASES_TMP}" "${JSON_OUT}" "${EPSILONS}" "${DELTA_NUMERS}" "${RATIO}" "${RUNS}" "${SIZE}" "${IDS}" "${HOST_INFO}" "${COMPILER_INFO}" "${CGAL_INFO}" "${GIT_SHA}" "${DATE_ISO}" <<'PY'
import json, sys, pathlib
tmp, out, epsilons, delta_numers, ratio, runs, size, ids, host, compiler, cgal, sha, date_iso = sys.argv[1:]
with open(tmp) as f:
    cases = json.load(f)
# Separate aggregates and per-id cases
aggs = [c for c in cases if c.get("_aggregate")]
per_id = [c for c in cases if not c.get("_aggregate")]
meta = {
    "date": date_iso,
    "host": host,
    "compiler": compiler,
    "cgal": cgal,
    "commit": sha,
    "epsilon_tiers": epsilons.split(),
    "delta_numers": delta_numers.split(),
    "runs": int(runs),
    "ratio": float(ratio),
    "size": size,
    "ids": ids.split(),
    "dataset": "derived via scripts/derive_benchmark_data.py (small ~100 pts IDs 11..20, large ~1000 pts IDs 21..30)",
    "method": "10-run Welford mean/stddev, Welch t 95% one-sided, SIMPLIFY_CORE_MS vs SQUISH_CORE_MS",
}
out_data = {
    "meta": meta,
    "aggregates": aggs,
    "cases": per_id,
}
with open(out, "w") as f:
    json.dump(out_data, f, indent=2)
print(f"Wrote {out} with {len(per_id)} cases and {len(aggs)} aggregates", file=sys.stderr)
PY

rm -f "${JSON_CASES_TMP}"

# Print markdown table
echo ""
echo "# SOTS vs SQUISH Benchmark"
echo ""
echo "Meta: host=\`${HOST_INFO}\` compiler=\`${COMPILER_INFO}\` CGAL=\`${CGAL_INFO}\` commit=${GIT_SHA} date=${DATE_ISO}"
echo ""
echo "Config: epsilons=(${EPSILONS}) delta_numers=(${DELTA_NUMERS}) size=${SIZE} ids=(${IDS}) runs=${RUNS} squish_ratio=${RATIO}"
echo ""
echo "| epsilon | delta | numer | SOTS mean ms ± std | SQUISH mean ms ± std | speedup (squish/SOTS) | Welch (95%) |"
echo "|---|---|---|---|---|---|---|"
# shellcheck disable=SC2059
printf "%s" "${MARKDOWN_ROWS}"
echo ""
echo "Speedup <1 means SQUISH faster (SOTS slower), >1 means SOTS faster."
echo ""
echo "CSV: ${CSV_OUT}"
echo "JSON: ${JSON_OUT}"
echo "TSV: ${TSV_OUT}"
echo ""
# Also print per-id detail from CSV for quick glance
echo "Per-ID details (first 10 rows):"
head -n 11 "${CSV_OUT}" | sed 's/^/| /;s/$/ |/;s/,/ | /g' | head -n 20

