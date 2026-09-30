| Algorithm | Mean Fréchet ↓ | Mean points ↓ | Mean core time ↓ | Speed vs SOTS |
|---|---|---|---|---|
| **SOTS** | 621 | **46.2 ✓** | 20.69 ms | — |
| DOTS | 477 | 75.0 | 0.24 ms | 87× faster |
| SQUISH | 454 | 89.8 | 0.046 ms | 450× faster |
| DP | 442 | 65.1 | 0.020 ms | 1,060× faster |

*Source: `results/fair-core/summary.csv` — 200 cases (20 trajectories × 5 ε × 2 bounds); values are case-balanced arithmetic means (10-run Welford). Bold ✓ = SOTS better (lower is better).*

# Simplification of Trajectory Streams

This repository contains the streaming delta-simplification algorithm from
[Simplification of Trajectory Streams](https://arxiv.org/abs/2503.23025), a
headless command-line program, a web visualizer, and comparison tooling
for trajectory-simplification baselines.

The project is research software. It is tested primarily on macOS arm64;
other platforms may work with equivalent CGAL, CMake, Julia, and C++
dependencies. Qt 6 Core is still required for the vendored DOTS baseline.

## Demo — Viewer walkthrough (1 min)

Short walkthrough of the web viewer on a real trajectory (Trajectory 1, 588 points, ε = 0.9, δ = 500): loading the stream, toggling original (gray) versus SOTS simplified result (green), scrubbing through streaming prefixes/steps, animating with Play/Pause and speed controls, and adjusting ε/δ. A second trajectory is loaded briefly to show scale.

https://github.com/user-attachments/assets/PLACEHOLDER

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
