# Simplification of Trajectory Streams

This repository contains the streaming delta-simplification algorithm from
[Simplification of Trajectory Streams](https://arxiv.org/abs/2503.23025), a
headless command-line program, an optional Qt viewer, and comparison tooling
for trajectory-simplification baselines.

The project is research software. It is tested primarily on macOS arm64;
other platforms may work with equivalent CGAL, Qt, CMake, Julia, and C++
dependencies.

## Repository layout

- `simplify.cpp`: headless implementation of the paper's algorithm.
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

- C++23 compiler (Clang 16+ or Apple Clang 15+)
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

On Ubuntu:

```bash
sudo apt update
sudo apt install build-essential cmake libcgal-dev libcgal-qt6-dev \
  qt6-base-dev julia python3-pip
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
include `-d DELTA`, `-e EPSILON`, `--dist`, and `--gui` on the GUI target.

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

The full benchmark is optional and can take hours. It requires the complete
raw dataset, Julia dependencies, and the `dots` target:

```bash
python3 scripts/benchmark.py --a 1 --b 1000 --workers 1 --resume
```

The benchmark writes generated files below `data/` and CSV output under
`results/`. Use `python3 scripts/benchmark.py --help` for resource limits,
timeouts, worker settings, and resume behavior.

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
