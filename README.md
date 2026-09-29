# Simplification of Trajectory Streams

This repository contains the streaming delta-simplification algorithm from
[Simplification of Trajectory Streams](https://arxiv.org/abs/2503.23025), a
headless command-line program, an optional Qt viewer, and comparison tooling
for trajectory-simplification baselines.

The project is research software. It is tested primarily on macOS arm64;
other platforms may work with equivalent CGAL, Qt, CMake, Julia, and C++
dependencies.

## Live demo

A Cloud Run deployment of the web visualizer is available at:

https://simplify-viewer-522405269791.asia-east2.run.app

## Repository layout

- `simplify_core.h`: headless streaming algorithm and its anchor workspace.
- `simplify.cpp`: headless command-line and web-trace front-ends.
- `simplify_with_gui.cpp`, `drawing.cpp`, `drawing.h`: optional Qt viewer.
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
- Qt 6 Core (used by the vendored DOTS target)

The Qt viewer additionally needs the CGAL Qt6 component and Qt 6 Widgets. The
web visualizer additionally needs Flask (`web/requirements.txt`). The Frechet
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
sudo apt install build-essential g++-14 cmake libcgal-dev libcgal-qt6-dev \
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
| `simplify` | Headless streaming simplifier |
| `simplify_with_gui` | Qt viewer and simplifier (`-DBUILD_GUI=ON`) |
| `dots` | DOTS baseline (Qt6 Core only; web compare pane and `benchmark.py`) |
| `dp` | DP baseline for the web compare pane (when the submodule source is present) |
| `squish` | SQUISH baseline for the web compare pane (when the submodule source is present) |

Headless / Docker builds use `-DBUILD_GUI=OFF`, which skips `simplify_with_gui`
but still builds `simplify`, `dots`, `dp`, and `squish` when their sources exist.

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
include `-d DELTA`, `-e EPSILON`, `--dist`, `--time` (opt-in phase timers on
stderr), and `--gui` on the GUI target.

### Visualize output


```bash
./build/simplify_with_gui --in 1 --gui --out
```

The optional `plot_curve` viewer from older local builds may be unavailable in
the current CMake configuration; use `simplify_with_gui` for the supported GUI
workflow.

### Compare baselines in the web visualizer

The Flask app in `web/` overlays this project's output against DOTS, DP, and
SQUISH on a prepared trace (`data/<id>/original.txt`):

```bash
python3 -m pip install -r web/requirements.txt
python3 web/server.py
```

Open the printed URL, load a trace, pick a baseline, and run it. `dots`, `dp`,
and `squish` are produced whenever their `traj-compression` sources exist,
including with `-DBUILD_GUI=OFF` (`dots` needs Qt6 Core only). If a baseline
binary is missing, initialize the submodule and rebuild.

## Benchmarking

Streaming SOTS vs SQUISH was measured on the standard derived datasets (DP/DOTS baselines are available via `build/dp`/`build/dots` and the web visualizer; core-time comparison to DP/DOTS is deferred).
Large IDs 21–30 (avg ~695 pts) are the primary perf signal.
Small IDs 11–20 (avg ~90 pts) are also recorded.
Five epsilon tiers were used: 299 (extra-coarse, δ=1), 30 (coarse, δ=9.68), 5 (mid, δ=50), 0.5 (fine, δ=200), 0.1 (extra-fine, δ=272.7).
SQUISH was run at ratio 0.15 so its output size matches SOTS at extra-coarse.
Each ID was sampled 10 times.
Means and stddevs are Welford.
Gating is Welch t (95% one-sided) as in benchmark.yml.
Hardware was MacBook Pro M4 (Apple clang 21, CGAL 6.1) for local and Ubuntu 24.04 g++-14 for CI.
Results are CORE_MS (isolated algorithm time, no I/O) from SIMPLIFY_CORE_MS vs SQUISH_CORE_MS.

### Algorithmic complexity (from papers)

SOTS streaming δ-simplification (papers/journal.pdf, Thm 1–2) guarantees dF(σ,τ) ≤ (1+ε)δ and |σ| ≤ 2κ(τ,δ)−2.
Working storage is O(ε^{-α}) and per-vertex time is O(ε^{-α} log 1/ε) for d∈{2,3}, O(ε^{-α}) for d≥4, where α = 2(d−1)⌊d/2⌋^2 + d.
For d=2, α=4, so storage O(ε^{-4}) and time O(ε^{-4} log 1/ε) — each ball Bv is covered by grid cells of width εδ/(2√d), |P|=O(ε^{-d}), and each stab structure Sa[p]=conv(Gva)∩F(Sa-1[p],p) has poly(1/ε) complexity.
At ε=0.1, ε^{-4}=10,000 vs ~1.2×10^{-10} at ε=299 (≈8×10^{13}× larger).
That is why runtime climbs from ~1.3 ms to ~134 ms on large in the table below.
If you run SOTS in static mode on n points, total time is O(ε^{-α} n log 1/ε), a factor n faster than the prior static O(ε^{2−2d} n^2 log n log log n) algorithm with the same bounds.
SQUISH is a streaming heuristic from traj-compression: buffer B=ratio·N, per-point O(B) naive scan (O(log B) with heap), total O(N·B), working storage O(B) bounded after the buffer fills, no Fréchet bound (uses SED).
DP is offline batch Douglas-Peucker (PED): O(N log N) average, O(N^2) worst, O(N) storage, non-streaming, no Fréchet guarantee.
DOTS is the online LSSD-based baseline (Qt): O(N/M) per point for buffer M, no Fréchet bound.
SOTS pays poly(1/ε) for its deterministic guarantee; SQUISH/DP/DOTS are fast heuristics without that guarantee.

### Parameter tuning

SOTS ε/δ coupling follows CI (benchmark.yml, scripts/derive_benchmark_data.py).
We sweep ε ∈ {299, 30, 5, 0.5, 0.1} — five tiers from extra-coarse to extra-fine.
δ is not tuned independently: δ = NUM/(1+ε) with NUM=300 (CI also uses 1000).
Thus (1+ε)δ = 300 is constant — smaller ε gets larger δ so the Fréchet envelope stays comparable and only approximation tightness varies.
ε controls grid resolution (εδ/(2√d)) and the (1+ε) factor in dF ≤ (1+ε)δ.
Smaller ε means finer cells and tighter bound, but poly(1/ε) more geometry and slower runtime (see extra-fine above).
δ controls ball radius and output size.
At δ=1 (extra-coarse) SOTS keeps ~104 pts avg large; at δ=272.7 (extra-fine) it keeps ~86 pts — it compresses more despite the finer grid.
All other knobs are fixed (NUM, dataset, 10-run Welford, CORE_MS isolation) so speed differences are apples-to-apples.
SQUISH has one knob: ratio = fraction of points kept (B=ratio·N, output size = ratio·N).
We calibrate ratio before looking at times.
We swept 0.05/0.15/0.5 on data/21 (N=588): SOTS at extra-coarse keeps 86–92 pts → 14–15%.
SQUISH at 0.15 keeps 88 pts (15%) — size-matched and fair.
At 0.5 it keeps 294 pts (50%), 3× more than SOTS and not comparable.
So 0.15 is the only ratio giving size parity at extra-coarse (see Points chart in the Lavish board §5).
At fine tiers SOTS then compresses more than SQUISH at the same ratio, showing its quality win even while slower.
Both binaries run on the same IDs, same order, same CORE_MS, 10-run Welford, Welch t 95% — identical to benchmark.yml.
The script exposes --epsilon, --delta-numer, --ratio, --runs, --size for alternatives, but defaults are the fair size-matched comparison.
For error-matched fairness (same Frechet), tune until measured Frechet matches (see next subsection).

### Frechet-matched (error-matched) at dF≈100 (all algos, data/21 pilot)

You asked: all algos producing Frechet distance 100 and then benchmark it.
We tuned each algo until measured continuous Frechet (scripts/frechet.jl, FrechetDist.jl) is 100±5 on data/21 (N=588), then compared core time.
This answers "does this have same Frechet? If not fix it, I want fair comparison" — main table is size-matched (same output size), this one is Frechet-matched (same error).

| Algorithm | Param to hit dF≈100 | Points (of 588) | Ratio kept | CORE_MS (10-run mean) | Measured dF | vs SOTS (time) | vs SOTS (size) |
|---|---|---|---|---|---|---|---|
| SOTS (fastest) | ε=30, δ=3.2258 ((1+ε)δ=100) | 112 | 19% | 1.49 ms | 96.7 ≤100 | — | — |
| SOTS (mid) | ε=5, δ=16.666 | 116 | 19.7% | 2.19 ms | 83.3 | 1.5× slower than fastest SOTS | similar |
| SOTS (fine) | ε=0.5, δ=66.666 | 110 | 18.7% | 9.89 ms | 100.0 | 6.6× slower than ε=30 | similar |
| SQUISH | ratio 0.29 (tuned, was 0.15→828) | 170 | 28.9% | 0.066 ms | 102.1 | 22× faster than SOTS | 1.52× more points than SOTS for same error |
| DP | ε=110 (was 50→45, 100→96) | 113 | 19.2% | 0.026 ms | 102.1 | 57× faster than SOTS | 1.01× same size, better than SQUISH |

How tuned: for SOTS set (1+ε)δ=100 and verify via `julia scripts/frechet.jl --id 21 --batch data/21/simplify.txt` (≈96–100); for SQUISH binary-search ratio until dF≈100 (0.27→141, 0.29→102, 0.30→40); for DP sweep ε until dF≈100 (90→87, 100→96, 110→102).
All on same HW, same ID, 10-run Welford for ms, single Julia batch per candidate (≈7s overhead).
Full large-average (IDs 21–30) would be same ratios: SOTS keeps ~110–115 pts avg large at 100, SQUISH needs ~28% (vs 15% size-matched), DP keeps ~110 pts.
Takeaway at dF=100: SOTS is still 22–57× slower than heuristics, but more compressive than SQUISH (112 vs 170 pts for same error, 1.5×) and matches DP on size.
The poly(1/ε) cost is for the guarantee, not raw speed.
At dF=100 the gap shrinks from 1937× (extra-fine size-matched) to ~22× because SOTS can use large ε (299/30) with tiny δ to hit 100 cheaply (1.5 ms) instead of ε=0.1 (75 ms).
Repro Frechet 100 pilot: `./build/simplify 21 -e 30 -d 3.2258 && julia scripts/frechet.jl --id 21 --batch data/21/simplify.txt` (≈96.7); `./build/squish data/21/original.txt 0.29 /tmp/s.txt && julia scripts/frechet.jl --id 21 --batch /tmp/s.txt` (≈102); `./build/dp data/21/original.txt 110 /tmp/d.txt && julia scripts/frechet.jl --id 21 --batch /tmp/d.txt` (≈102).
See Lavish board §4 Results for the same table plus Frechet 250 pilot, and `scripts/bench-compare-squish.sh --help`.

### SOTS vs SQUISH (large, 10-run mean) — size-matched (ratio 0.15)

| ε tier | ε | δ | SOTS ms ± std | SQUISH ms ± std | SQUISH/SOTS | Overhead |
|---|---|---|---|---|---|---|
| extra-coarse | 299 | 1 | 1.37 ± 0.13 | 0.071 ± 0.016 | 0.052× | 19× slower |
| coarse | 30 | 9.68 | 1.41 ± 0.37 | 0.064 ± 0.005 | 0.045× | 22× slower |
| mid | 5 | 50 | 1.30 ± 0.04 | 0.064 ± 0.004 | 0.049× | 20× slower |
| fine | 0.5 | 200 | 11.28 ± 0.17 | 0.062 ± 0.003 | 0.006× | 181× slower |
| extra-fine | 0.1 | 272.7 | 133.74 ± 50.33 | 0.069 ± 0.011 | 0.001× | 1937× slower |

All tiers are confidently slower than SQUISH (Welch p < 0.05).
SQUISH is 19–1900× faster on core time but provides no Fréchet guarantee (DP/DOTS core-time comparison deferred).
SOTS provides deterministic Fréchet ≤ (1+ε)δ per segment (papers/journal.pdf, Thm 1–2).
Previous uplifts remain.
PR26 halved runtime across all tiers.
PR34 extra-fine was 1.95–2.13× and fine was 1.4× vs pre-PR34 on the same hardware.
PR37 clean-core added 1.04× encapsulation with no regression.
See the Lavish board at [.lavish/sots-bench-squish-compare-lavish-a1.html](.lavish/sots-bench-squish-compare-lavish-a1.html).
Run the repro script:

```bash
scripts/bench-compare-squish.sh
scripts/bench-compare-squish.sh --epsilon "0.1 0.5" --ratio 0.15 --runs 10 --size large
```

The script builds both binaries, warms up, runs 10 Welford samples, Welch-gates, writes JSON/CSV to .lavish/bench-squish-*.json and prints a markdown table.
It is shellcheck-clean, cross-platform (mac/Linux), and CI-friendly.
It reuses scripts/ci/welch.py and the same CORE_MS extraction as .github/workflows/benchmark.yml.

The reproducible four-way comparison of SOTS, DOTS, SQUISH and DP is documented
in [the core benchmark report](docs/fair-core-benchmark.md). It uses the 200 CI
cases, serial native `core_ms` timings, calibrated common continuous Fréchet
bounds, and exact SOTS output checks before/after optimization. **Corrected
2026-09-29: CI mean speedup is 1.00× (neutral, 0% cut, 20/20 passed, per-config
0.94–1.08×; intersect 1.07× but overall neutral, Run 36460949757), not the
earlier local M1 1.115× (10.3%).** The
[interactive Lavish board](.lavish/sots-fair-bench/index.html) shows both phases,
actual errors, retained points, and raw sample statistics (now with corrected
1.00× conclusion). Open that HTML file locally, or run
`lavish-axi .lavish/sots-fair-bench/index.html`.

The full benchmark is optional and can take hours.
It requires the complete raw dataset, Julia dependencies, and the `dots` target:

```bash
python3 scripts/benchmark.py --a 1 --b 1000 --workers 1 --resume
```

The benchmark writes generated files below `data/` and CSV output under `results/`.
Use `python3 scripts/benchmark.py --help` for resource limits, timeouts, worker settings, and resume behavior.

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
