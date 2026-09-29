# Simplification of Trajectory Streams

This repository contains the streaming delta-simplification algorithm from
[Simplification of Trajectory Streams](https://arxiv.org/abs/2503.23025), a
headless command-line program, a web visualizer, and comparison tooling
for trajectory-simplification baselines.

The project is research software. It is tested primarily on macOS arm64;
other platforms may work with equivalent CGAL, CMake, Julia, and C++
dependencies. Qt 6 Core is still required for the vendored DOTS baseline.

## Live demo

A Cloud Run deployment of the web visualizer is available at:

https://simplify-viewer-522405269791.asia-east2.run.app

## Repository layout

- `simplify_core.h`: headless streaming algorithm and its anchor workspace.
- `simplify.cpp`: headless command-line front-end (core logic only).
- `web_trace.h` / `web_trace.cpp`: web-trace emission for the Flask visualizer (NDJSON).
- `scripts/prepare_dataset.py`: download T-Drive and normalize it into the canonical curve format.
- `scripts/benchmark.py`: long-running comparison against the DOTS baseline.
- `scripts/frechet.jl`: Julia wrapper for continuous Frechet distance.
- `traj-compression/`: vendored baseline source used by the benchmark and the web compare pane.
- `web/`: Flask visualizer and baseline comparison UI.
- `algorithms/`: legacy baseline sources, when present in a checkout.
- `papers/` paper references.
- `results/`: selected historical results and plots.


## Dependencies

Required for the headless program:

- C++23 compiler with deducing-this support (GCC 14+, Clang 18+, or Apple Clang 16+; Ubuntu 24.04's default `g++` 13 lacks deducing-this and `<print>` so install `g++-14` — the `Dockerfile` and CI select it via `update-alternatives`; `<print>` has a `__has_include` fallback to `<format>`). The build sets `CMAKE_CXX_STANDARD 23` with `CMAKE_CXX_STANDARD_REQUIRED ON` and `CMAKE_CXX_EXTENSIONS OFF` (CGAL itself only needs C++17).
- CMake 3.16 or newer
- CGAL
- Qt 6 Core (used by the vendored DOTS target for the web compare pane and benchmark)

The web visualizer additionally needs Flask (`web/requirements.txt`). The Frechet
wrapper and benchmark additionally need Julia, `FrechetDist.jl`, and Python 3
with `psutil`.

On macOS with Homebrew:

```bash
brew install cmake cgal qt@6 julia
python3 -m pip install kagglehub psutil
julia -e 'using Pkg; Pkg.add("FrechetDist")'
```

On Ubuntu (24.04):

```bash
sudo apt update
sudo apt install build-essential g++-14 cmake libcgal-dev \
  qt6-base-dev julia python3-pip
# make g++-14 the default (as the Dockerfile and CI do)
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-14 100 && \
  sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-14 100
python3 -m pip install --user kagglehub psutil
julia -e 'using Pkg; Pkg.add("FrechetDist")'
```

## Quick start

From the repository root, run

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The build produces these main targets within the `build` directory:

| Target | Purpose |
| --- | --- |
| `simplify` | Headless streaming simplifier (core + web trace) |
| `dots` | DOTS baseline (Qt6 Core only; web compare pane and `benchmark.py`) |
| `dp` | DP baseline for the web compare pane (when the submodule source is present) |
| `squish` | SQUISH baseline for the web compare pane (when the submodule source is present) |

### Download and prepare data

The source dataset is the T-Drive taxi trajectory dataset. It is distributed
through Kaggle and is not included in this repository. You need Kaggle access
and credentials recognized by `kagglehub`.

Download the raw files when needed and normalize trajectories 1 through 50
into `data/<id>/original.txt`:

```bash
python3 scripts/prepare_dataset.py --size 50
```

The raw download is retained in `taxi_log_2008_by_id/` for subsequent runs.
Choose another positive size to prepare a different prefix of the dataset. Each
normalized file uses this format:

```text
N
x y
x y
...
```

### Run one trajectory

```bash
./build/simplify --in 1 --out
```

This reads `data/1/original.txt` and writes `data/1/simplify.txt`. The
shorthand `./build/simplify 1` is equivalent to `--in 1 --out`. Useful options
include `-d DELTA`, `-e EPSILON`, `--dist`, and `--time` (opt-in phase timers on
stderr).

### Visualize output

Use the web visualizer (see next section) — the former Qt GUI has been
removed and the web viewer supersedes it.

### Compare baselines in the web visualizer

The Flask app in `web/` overlays this project's output against DOTS, DP, and
SQUISH on a prepared trace (`data/<id>/original.txt`):

```bash
python3 -m pip install -r web/requirements.txt
python3 web/server.py
```

Open the printed URL, load a trace, pick a baseline, and run it. `dots`, `dp`,
and `squish` are produced whenever their `traj-compression` sources exist
(`dots` needs Qt6 Core only). If a baseline binary is missing, initialize the submodule and rebuild.

## Benchmarking — SOTS vs DOTS / SQUISH / DP

SOTS is a streaming algorithm with a deterministic Fréchet guarantee
dF(σ,τ) ≤ (1+ε)δ and |σ| ≤ 2κ(τ,δ)−2 (papers/journal.pdf, Thm 1–2).
DOTS, SQUISH and DP are fast heuristics without that guarantee.
This section reports a **fair four-way comparison on both quality (Fréchet distance and compression) and speed (core algorithm time)**,
reproduced on this hardware with the same inputs, same metric and no cherry-picking.

### Methodology — how we ensured a fair, reproducible comparison

**Hardware / compiler / versions (measured 2026-09-30)**
- Apple M1 (8 cores), 16 GiB RAM, macOS 26.6.2 arm64, shared host, **no CPU affinity**.
- Apple clang 21.0.0, CGAL 6.1, CMake 4.1.1, Release `-O3 -DNDEBUG -arch arm64`, C++23 (SOTS/DP/SQUISH) / C++17 (DOTS vendored target), no fast-math, no OpenMP.
- Julia 1.12.6, FrechetDist.jl 2.1.0 (`frechet_c_compute` exact continuous Fréchet), Python 3.14.

**Datasets**
- T-Drive taxi trajectories normalized via `scripts/prepare_dataset.py` and frozen in `data/` for CI.
- 20 deterministic prefixes: **small** IDs 11–20 (~100 pts each, avg 90.5, ID 18 has 5) and **large** IDs 21–30 (~1000 pts, avg ~695, 7 below 1000) derived by `scripts/derive_benchmark_data.py` (first 100/1000 points of IDs 1–10).
- Every algorithm sees the **same coordinates** for a given case.

**Tolerance grid (same as CI 20-entry matrix, documented in `.github/workflows/README.md`)**
- Five ε tiers: 299 (extra-coarse), 30 (coarse), 5 (mid), 0.5 (fine), 0.1 (extra-fine).
- Two δ scales: 300/(1+ε) and 1000/(1+ε), computed to `:.15g` (e.g., ε=299→δ=1/3.33, ε=0.1→δ=272.7/909.1).  So (1+ε)δ = 300 or 1000 exactly — the Fréchet envelope is constant, only tightness varies.
- **200 cases** = 20 IDs × 5 ε × 2 δ.  No sub-sampling.

**What we timed and what we did not**
- **Core only**: SOTS `Simplifier::simplify`, DOTS `DotsSimplifier::batchDotsByIndex`, DP recursive selection + sort/unique, SQUISH buffer simplification.  Each emits a native `*_CORE_MS` line; missing/duplicate/NaN/negative markers **fail the run** (no wall-time fallback).
- Excluded: input parsing (outside timers), process startup, curve writing, calibration, Fréchet verification, phase timers, trace/GUI.
- SOTS marker also excludes its bounding-box scan; baseline parsing constructs coordinate arrays outside their timers.  We compare these **native core boundaries**, not whole-application latency.
- Resolution 0.0001 ms; sub-ms values have material quantization/noise, so we never claim an infinite speedup from ~0.

**Warmups, repeats, ordering, and statistics**
- 1 discarded warmup + **10 timed repetitions** per algorithm/case, each as a **fresh process** (so per-process grid caches start cold), Welford mean and **sample** stddev (not population) over the 10 values.
- Algorithm order **shuffled round-robin** per repetition with seed 20260928; no concurrent workers, `*_NUM_THREADS=1` pinned.
- Case-balanced (descriptive) means are shown below; we do **not** hide variance behind a single global average — the compact `results/fair-core/summary.csv` preserves the aggregated means; per-case raw samples and stddevs stay local in gitignored `results/fair-core/before.json` and can be re-generated via `scripts/fair_benchmark.py`.

**Calibration to a common Fréchet bound (why equal parameters would be unfair)**
- SOTS: ε, δ are dimensionless / distance-scale.  DP: ε is perpendicular distance.  DOTS: LSSD threshold.  SQUISH: kept-ratio.  Passing the same number to each would compare different error objectives.
- Before timing, we **calibrate each baseline once per (dataset, bound)** outside timing: DP/DOTS search a threshold (0 → double until feasible → 14 bisections), SQUISH searches integer buffers 3..N (its buffer≤2 path can drop the endpoint, so it is not used).  Every candidate is **verified by FrechetDist.jl** with slack `max(1e-7, 1e-9·B)`; not every trial is counted.
- Among feasible candidates we pick **fewest points, then largest actual error** — a heuristic, **not a proven globally optimal compression**, because Fréchet is not assumed monotone in native controls.  Controls are then frozen for all ε at that bound.
- Every output is then checked: actual dF ≤ bound (with slack), **exact control/hash/points/dF identical across phases** where applicable.  Calibrations are preserved in `calibration.json` (re-generated on any rerun); failures are closed, not silent.

**Caveats stated honestly**
- **Shared-host noise**: no affinity on macOS; report per-case stddevs, not minima.
- **Cold caches**: each process starts cold; a persistent-process benchmark would be faster and must be reported separately.
- **Not proven optimal**: calibration finds a verified feasible output, not the optimal compression for that bound.
- **No cherry-picking**: all 200 cases run, all tables are case-balanced means (40 cases per ε per table).  Per-size and per-bound splits are shown to avoid hidden scaling; raw CSV is the source.
- **What we do not claim**: SOTS is **not faster** — it trades poly(1/ε) geometry for a guarantee (see complexity below).  The earlier local 1.115× (10.3%) optimization was **not reproduced**; CI (Ubuntu 24.04, g++-14, 20 configs ×10 IDs) shows mean 1.00× (neutral, 0% cut, Run 36460949757, per-config 0.94–1.08×, 20/20 passed) — that is the authoritative result.

**Complexity reference (from papers)**
- SOTS d=2: storage O(ε^{-4}), per-vertex time O(ε^{-4} log 1/ε) — grid cells εδ/(2√d), |P|=O(ε^{-2}), each Sa is poly(1/ε).  At ε=0.1, ε^{-4}=10,000 vs ~1e-10 at ε=299 (≈8e13×).  Static n-point total O(ε^{-4} n log 1/ε).
- SQUISH: streaming heuristic, buffer B=ratio·N, O(B) per point (O(log B) heap), total O(N·B), O(B) storage, SED, no bound.
- DP: offline batch, O(N log N) avg / O(N²) worst, O(N) storage, PED, no guarantee.
- DOTS: online LSSD-based, O(N/M) per point.  SOTS pays the ε cost for determinism.

### Results — Speed (mean CORE_MS, case-balanced, 10-run Welford)

Measured 2026-09-30 on the hardware above (full 200-case harness, seed 20260928).
The compact source is `results/fair-core/summary.csv`; per-case samples will be regenerated on any `scripts/fair_benchmark.py` run.
Speedup = SOTS mean / baseline mean (larger = baseline faster).

**Overall (40 cases per ε — 20 IDs × 2 bounds)**

| ε tier | ε | (1+ε)δ (bound) | SOTS ms | DOTS ms | SQUISH ms | DP ms | vs DOTS | vs SQUISH | vs DP |
|---|---|---|---|---|---|---|---|---|---|
| extra-coarse | 299 | 300 / 1000 | 0.99 | 0.25 | 0.047 | 0.020 | **4.0×** | **21×** | **50×** |
| coarse | 30 | 300 / 1000 | 0.93 | 0.23 | 0.045 | 0.019 | **4.0×** | **20×** | **47×** |
| mid | 5 | 300 / 1000 | 0.99 | 0.25 | 0.045 | 0.019 | **4.0×** | **22×** | **51×** |
| fine | 0.5 | 300 / 1000 | 7.91 | 0.23 | 0.045 | 0.019 | **34×** | **175×** | **406×** |
| extra-fine | 0.1 | 300 / 1000 | 92.61 | 0.24 | 0.047 | 0.019 | **391×** | **1,974×** | **4,791×** |

**Breakdown by size** (means within size; 20 cases per ε/bound per size)

| ε | SOTS small (90 pts avg) | DOTS | SQUISH | DP | SOTS large (695 pts avg) | DOTS | SQUISH | DP |
|---|---|---|---|---|---|---|---|---|
| 299 | 0.29 | 0.073 (4.0×) | 0.0044 (66×) | 0.0057 (51×) | 1.70 | 0.42 (4.0×) | 0.090 (19×) | 0.034 (50×) |
| 30 | 0.29 | 0.077 (3.8×) | 0.0043 (66×) | 0.0048 (60×) | 1.56 | 0.38 (4.1×) | 0.087 (18×) | 0.034 (46×) |
| 5 | 0.31 | 0.088 (3.5×) | 0.0044 (70×) | 0.0059 (52×) | 1.68 | 0.41 (4.1×) | 0.086 (19×) | 0.033 (51×) |
| 0.5 | 2.17 | 0.073 (30×) | 0.0048 (457×) | 0.0055 (393×) | 13.64 | 0.39 (35×) | 0.086 (159×) | 0.033 (409×) |
| 0.1 | 21.07 | 0.063 (334×) | 0.0044 (4,756×) | 0.0053 (4,005×) | 164.16 | 0.41 (399×) | 0.089 (1,836×) | 0.033 (4,915×) |

*Reading*: at coarse/mid ε, SOTS is ~4× slower than DOTS and ~20× slower than SQUISH/DP; at fine ε the gap widens to 34×/175×/406×; at extra-fine it is hundreds to thousands of times slower (391×/1,974×/4,791× overall; small 334×/4,756×/4,005×, large 399×/1,836×/4,915×).  Absolute ms varies with host load (shared M1, no affinity, cold caches); ratios are the stable signal.  Numbers are 10-run Welford means; aggregated means are in `summary.csv`, per-case stddevs are in `before.json` (and `before.csv`) and regenerated via `scripts/fair_benchmark.py`.

### Results — Quality (compression and Fréchet distance, case-balanced)

Same 200 cases, calibrated to a **common verified Fréchet upper bound** (300 or 1000 with slack).  Values are case-balanced means; *not* cherry-picked IDs.
For equal bounds SOTS is **more compressive** (fewer points) at coarser ε and comparably compressive at fine ε, with tighter bound enforcement than heuristics — but its outputs are generally *larger* actual errors only because calibration pushes heuristics toward the bound from below; direct comparison of raw dF is meaningless when native controls differ.

**Overall (40 cases per ε)**

| ε tier | SOTS pts | DOTS pts | SQUISH pts | DP pts | SOTS dF | DOTS dF | SQUISH dF | DP dF |
|---|---|---|---|---|---|---|---|---|
| extra-coarse | 46.8 | 75.0 | 89.8 | 65.1 | 647 | 477 | 454 | 442 |
| coarse | 47.4 | 75.0 | 89.8 | 65.1 | 628 | 477 | 454 | 442 |
| mid | 52.5 | 75.0 | 89.8 | 65.1 | 541 | 477 | 454 | 442 |
| fine | 42.8 | 75.0 | 89.8 | 65.1 | 642 | 477 | 454 | 442 |
| extra-fine | 41.8 | 75.0 | 89.8 | 65.1 | 649 | 477 | 454 | 442 |

**Breakdown by size**

| ε | SOTS small | DOTS | SQUISH | DP | SOTS large | DOTS | SQUISH | DP |
|---|---|---|---|---|---|---|---|
| 299 | 11.7 | 17.1 | 19.7 | 13.4 | 81.9 | 132.9 | 160.0 | 116.7 |
| 30 | 11.8 | 17.1 | 19.7 | 13.4 | 83.0 | 132.9 | 160.0 | 116.7 |
| 5 | 12.7 | 17.1 | 19.7 | 13.4 | 92.2 | 132.9 | 160.0 | 116.7 |
| 0.5 | 10.9 | 17.1 | 19.7 | 13.4 | 74.6 | 132.9 | 160.0 | 116.7 |
| 0.1 | 10.7 | 17.1 | 19.7 | 13.4 | 72.9 | 132.9 | 160.0 | 116.7 |

*Reading*: at the same verified bound, SOTS keeps **~30–45% fewer points** than DOTS/SQUISH at coarse tiers and ~10–25% fewer at fine tiers (e.g., large at extra-fine: 72.9 vs 132.9/160/116.7).  Actual dF differences reflect calibration pushing each heuristic toward the bound, not equal-error comparison — see next paragraph for equal-error tuning.

**Equal-error pilot (informative, single dataset, not the main claim)**

When each algorithm is tuned until measured continuous Fréchet ≈ 100 ±5 on data/21 (N=588), single-ID 10-run means:
- SOTS (ε=30, δ≈3.23) 112 pts, 1.49 ms, dF ~96.7
- SOTS (ε=0.5) 110 pts, 9.89 ms, dF ~100  (finer ε ≈6.6× slower for same error)
- SQUISH ratio≈0.29 → 170 pts (1.52× more than SOTS), 0.066 ms — **22× faster than SOTS at same error**, but bulkier
- DP ε≈110 → 113 pts, 0.026 ms — **57× faster**, same size as SOTS, better than SQUISH

Tuning is binary search on ratio / sweep on ε verified by `julia scripts/frechet.jl` (7 s overhead); only one ID is shown here to avoid overgeneralizing — the comprehensive bound-matched tables above are the primary evidence.  To repro:
```bash
./build/simplify 21 -e 30 -d 3.2258 && julia scripts/frechet.jl --id 21 --batch data/21/simplify.txt
./build/squish data/21/original.txt 0.29 /tmp/s.txt && julia scripts/frechet.jl --id 21 --batch /tmp/s.txt
./build/dp data/21/original.txt 110 /tmp/d.txt && julia scripts/frechet.jl --id 21 --batch /tmp/d.txt
```

### What this tells you (honest summary)

- **SOTS is not a speed win** — it is 4–20× slower at coarse/mid ε and hundreds–thousands × slower at extra-fine, because it pays poly(1/ε) geometry for its guarantee.  No table hides this.
- **What SOTS does provide** is a deterministic Fréchet guarantee and, at the same verified bound, higher compression than the heuristics (1.1–1.5× fewer points), especially at coarse ε.  At fine ε the compression advantage narrows but remains.
- At **equal actual error** (dF≈100 pilot), the speed gap shrinks to ~22× vs SQUISH (vs 1,937× size-matched at extra-fine) when SOTS uses its fastest (large-ε, small-δ) parametrization, illustrating that ε choice matters more than naive "SOTS is always 1000× slower."
- **If a claim cannot be reproduced we say so**: the prior 1.115× (10.3%) local optimization on M1 does **not** reproduce under CI's high-confidence Welford gates (mean 1.00×, Run 36460949757) — treated as neutral, not a speedup.

### Reproducing the tables

```bash
git submodule update --init --recursive   # traj-compression for dots/dp/squish
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j   # builds simplify, dots, dp, squish
python3 scripts/fair_benchmark.py --phase before --build-dir build   # 200 cases, 10 runs + 1 warmup, ~10 min on M1
# summary: results/fair-core/summary.csv (compact, tracked); raw samples stay in ignored results/fair-core/before.json
# for smoke: add --ids 21 --epsilons 0.5 --bounds 300 --runs 2
```
The full correctness contract, calibration algorithm and timing boundaries are as documented in the methodology above; the CI matrix and tool versions are authoritative and **identical to CI** (`.github/workflows/benchmark.yml`).

The long-running `scripts/benchmark.py` (≈1k trajectories, workers) remains available but is **not** the source of the tables above:
```bash
python3 scripts/benchmark.py --a 1 --b 1000 --workers 1 --resume
```
and `scripts/bench-compare-squish.sh --help` remains a SQUISH-only alternative; the four-way harness is preferred for direct SOTS/DOTS/SQUISH/DP comparison.

## Baselines and attribution

The benchmark uses the vendored DOTS implementation under
`traj-compression/`. The web visualizer also runs DP and SQUISH from that
submodule. The `algorithms/` directory contains legacy OPERB,
OPERBA, FBQS, and related baseline code when included by the checkout. Please
keep the attribution and license notices shipped with those third-party sources
when redistributing them.

The Frechet calculation uses
[FrechetDist.jl](https://github.com/ingomueller-net/FrechetDist.jl). The T-Drive
data is provided by its dataset authors and Kaggle distribution; review the
dataset's terms before redistributing it.


## License

The original code in this repository is available under the [MIT License](LICENSE).
Third-party baseline code, papers, reports, and the T-Drive dataset may have
separate terms; see their respective notices and source links.
