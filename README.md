# Simplification of Trajectory Streams

This repo implements the delta-simplification from the paper
[Simplification of Trajectory Streams](https://arxiv.org/abs/2503.23025). It
includes a Qt viewer and a CLI to simplify trajectories and optionally
compute (approximate) Fréchet distance between the original and simplified
polylines. It also benchmarks this streaming algorithm against several
classical baselines (DP, DOTS, SQUISH).

Tested primarily on macOS (arm64). Other platforms may work with the right
dependencies.

Demo video of this algorithm:
<https://github.com/user-attachments/assets/7593b5d8-ca68-4b7b-af6f-de76cd1c7689>

Visual comparison between different algorithms (the algorithm implmenete is in red):
<img width="724" height="725" alt="Comparison screenshot" src="https://github.com/user-attachments/assets/d5bc2041-b0b6-47a6-a261-7ef5a7eac757" />

---

## Prerequisites

- C++ toolchain with **C++23** support (Clang 16+ or Apple Clang 15+)
- CMake ≥ 3.16
- **CGAL** (with Qt6 support, for the viewer)
- **Qt 6 Widgets** (QtCore, QtGui, QtWidgets)
- **Julia** with the `FrechetDist` package installed (required for the
  `./scripts/frechet` wrapper used by `--dist` and the benchmark)

On macOS with Homebrew:

```bash
brew install cmake cgal qt@6 julia
```

For Julia:

```bash
julia -e 'using Pkg; Pkg.add("FrechetDist")'
```

---

## What's in the repo

```
.
├── CMakeLists.txt
├── README.md
│
├── simplify.cpp               # Headless main algorithm (this paper)
├── simplify_with_gui.cpp      # Main algorithm with the Qt viewer
├── simplify_old.cpp           # Older snapshot (kept for reference)
├── drawing.cpp / drawing.h    # Shared Qt viewer widget
│
├── scripts/
│   ├── normalize.cpp          # Convert raw taxi CSVs to data/<id>/original.txt
│   ├── plot_curve.cpp         # Multi-curve overlay viewer for data/<id>/
│   ├── benchmark.py           # Overnight benchmark of simplify vs DP/DOTS/SQUISH
│   ├── download_dataset.py    # Helper to fetch the T-Drive dataset from Kaggle
│   ├── clean_data.sh          # Strip non-original.txt files from data/<id>/
│   └── frechet                # Julia wrapper around FrechetDist.jl (used by
│                              #   the C++ --dist flag and benchmark.py)
│
├── algorithms/                # Subtree with OPERB/OPERBA/FBQS sources (vendored from
│   └── ...                    #   Trajectory-Simplification-Algorithm, used by algorithms/main)
│
├── traj-compression/          # Subtree of baseline algorithm sources
│   ├── batch/DP/              #   - DP_adapted (Douglas–Peucker)
│   ├── online/DOTS/           #   - DOTS_adapted
│   └── online/SQUISH/         #   - SQUISH_adapted
│
├── data/                      # Sample trajectories (102 ids committed; rest is generated)
│   ├── <id>/                  # Per-id input + outputs:
│   │   ├── original.txt       #   - N points "x y" each, after normalization
│   │   ├── simplify.txt      #   - the streaming algorithm's output
│   │   ├── DOTS.txt           #   - DOTS output (if benchmarked)
│   │   ├── DP.txt             #   - DP output
│   │   ├── SQUISH.txt         #   - SQUISH output
│   │   └── simplify_against_<baseline>_<e>_<d>.txt   # simplify tuned per baseline
│   ├── taxi/                  # Normalized raw inputs (one .txt per id)
│   └── taxi_simplified/       # Where the legacy algorithms/main wrote its outputs
│
├── report/                    # UROP progress report and figures
└── icons/, resources/         # Qt resources used by the viewer
```

---

## Build

```bash
mkdir -p release
cd release
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . -j
```

Targets produced inside `release/`:

| Target | Purpose |
|---|---|
| `simplify` | Headless main algorithm (this paper). |
| `simplify_with_gui` | Main algorithm with the Qt viewer. |
| `plot_curve` | Overlay viewer for every curve in `data/<id>/`. |
| `normalize` | Convert raw T-Drive CSVs to normalized `data/<id>/original.txt`. |
| `DP_adapted` | Douglas-Peucker baseline (vendored from `traj-compression/batch/DP/`). |
| `DOTS_adapted` | DOTS baseline (vendored from `traj-compression/online/DOTS/`). |
| `SQUISH_adapted` | SQUISH baseline (vendored from `traj-compression/online/SQUISH/`). |

The baseline binaries used by `benchmark.py` (`DP_adapted`, `DOTS_adapted`,
`SQUISH_adapted`) are built by the same `cmake --build .` invocation — no
extra per-subdir `make` step is needed. They live next to `simplify` so
the benchmark can find them in one place.

The `algorithms/` subdirectory is a standalone C++ project (its own
`Makefile`); build it separately if you want OPERB/OPERBA/FBQS as baselines
(see [Baseline algorithms](#baseline-algorithms)).

---

## Complete workflow: from a fresh clone to a working program

These steps reproduce everything end-to-end. Steps 1–3 are required; the
rest is the "full" benchmark path.

### 1. Install system and language dependencies

```bash
# macOS
brew install cmake cgal qt@6 julia
julia -e 'using Pkg; Pkg.add("FrechetDist")'
```

### 2. Clone and build the C++ scripts

```bash
git clone <this-repo>
cd Simplification-of-Trajectory-Streams
mkdir -p release
cd release
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . -j
cd ..
```

This produces `release/simplify`, `release/plot_curve`, `release/normalize`, etc.

### 3. Make sure the `data/` directory is ready

A partial set of 102 trajectories is already committed under `data/`. Each
`<id>/` directory contains an `original.txt` (the input) and zero or more
`*.txt` files (the outputs). You can use these directly.

If you want the **full** T-Drive dataset, fetch it with the helper script:

```bash
pip3 install kagglehub
python3 download_dataset.py      # downloads to ~/.cache/kagglehub/...
```

Then run the normalizer to produce `data/<id>/original.txt` for each
trajectory:

```bash
./release/normalize --all
# or a single id:
./release/normalize 16
```

The normalizer reads the last two comma-separated fields of each line in
`taxi_log_2008_by_id/<id>.txt`, scales the trajectory to roughly
`[-10000, 10000]` (preserving aspect ratio), and writes
`data/<id>/original.txt` in the format `N\nx y\nx y\n...`.

### 4. Run a single simplification

```bash
# Headless: read original.txt, write data/<id>/simplify.txt
./release/simplify --in 1 --out

# With the Qt viewer
./release/simplify_with_gui --in 1 --gui

# GUI + write output + compute Fréchet distance after the viewer closes
./release/simplify_with_gui --in 1 --gui --out --dist

# Override delta/epsilon
./release/simplify_with_gui --in 1 -d 1500 -e 0.5 --gui
```

Shorthand: `simplify <id> [flags]` is equivalent to
`--in <id> --out [flags]`.

**Path resolution.** The `simplify_*` binaries locate the repository root by
walking up from `argv[0]` until they find a `data/` directory (up to 5
levels), then fall back to the current working directory. So you can run
them from anywhere — no need to `cd release` first.

**Output layout** for `data/<id>/`:

```
data/<id>/
  original.txt          # input (created by scripts/normalize)
  simplify.txt          # streaming algorithm output (after --out)
```

### 5. Visualize a result

```bash
./release/plot_curve 1           # shows original + every *.txt in data/1/
```

Click legend entries to toggle curves. The viewer auto-loads anything that
matches the canonical `N x y` format.

### 6. Compute Fréchet distance

The `--dist` flag on `simplify` shells out to the Julia wrapper `frechet`
(path: `scripts/frechet`). You can also call it directly:

```bash
./scripts/frechet <id>

# Or with an explicit simplified path:
./scripts/frechet --in 1 --path data/1/simplify.txt
```

### 7. Run the full benchmark

`scripts/benchmark.py` is the long-running comparison. It runs DP, DOTS, and
SQUISH on each id, then sweeps `simplify` over a few
`(epsilon, delta)` values tuned per baseline, and records point counts and
Fréchet distances into `compare_points.csv`.

The baseline binaries (`DP_adapted`, `DOTS_adapted`, `SQUISH_adapted`) are
built into `release/` by the top-level CMake build above, so the benchmark
script does not need any extra build step.

```bash
python3 scripts/benchmark.py            # all ids
python3 scripts/benchmark.py --a 1 --b 3 --resume   # pick ids, resume from CSV
```

The benchmark has built-in safety nets: a 1.5 GB RSS cap per invocation
and per-stage timeouts. A run that exceeds either is recorded as `-1` in
the CSV and the loop continues.

### 8. Clean up output files (optional)

`clean_data.sh` keeps only `original.txt` in each `data/<id>/`, removing
all generated `simplify_*.txt`, `DP.txt`, `DOTS.txt`, etc. Useful before
re-running a benchmark.

```bash
./clean_data.sh            # interactive prompt
./clean_data.sh --force    # delete immediately
./clean_data.sh --dry-run  # show what would be deleted
```

---

## CLI reference: `simplify`

```
simplify [options] [<id>]

Options:
  --in <id>        Read input from data/<id>/original.txt
  --out            Write output to data/<id>/simplify.txt
  --dist           After output (or after the viewer closes), run frechet
                   to print Fréchet distance between original and simplified
  --gui            Show the Qt viewer while simplifying
  --labels         Show vertex index labels in the viewer
  -d <delta>       Override DELTA (default 200)
  -e <epsilon>     Override EPSILON (default 0.6)
  -F / -G / -S     Debug polygon overlays in the viewer
                   (F, G, S refer to the structures in the paper)
  -h, --help       Print this help and exit

Shorthand:
  simplify <id> [flags]   ==  --in <id> --out [flags]
```

---

## CLI reference: `scripts/frechet`

```
scripts/frechet [<id>] [--in <id>] [--path <file>] [--raw]

  <id>             Trajectory id (uses data/<id>/original.txt vs
                   data/<id>/simplify.txt)
  --in <id>        Same as positional id
  --path <file>    Override the simplified curve to compare against
  --raw            Print only the raw distance number
```

The wrapper is a Julia script that calls `frechet_c_compute` from
`FrechetDist.jl` internally.

---

## Baseline algorithms

The streaming algorithm is benchmarked against three groups of baselines:

| Group | Where | How to build | Used by |
|---|---|---|---|
| DP, DOTS, SQUISH | `traj-compression/` sources, `release/` binaries | `cmake --build release/` | `benchmark.py` |
| OPERB, OPERBA, FBQS | `algorithms/` | `cd algorithms && g++ -O2 -std=c++17 -Iinc main.cpp -o main && ./main` | manual / `plot_curve` |

The `algorithms/` subtree is vendored from
[Trajectory-Simplification-Algorithm](https://github.com/MingjiHan99/Trajectory-Simplification-Algorithm)
and modified for our I/O format. The legacy runner is invoked as:

```bash
g++ -O2 -std=c++17 -Iinc algorithms/main.cpp -o algorithms/main
./algorithms/main <id> [error_bound]
# Example: ./algorithms/main 1 100
```

This reads `data/<id>/original.txt`, runs all four algorithms, and writes
`data/<id>/dp_simplified.txt`, `operb_simplified.txt`, `operba_simplified.txt`,
`fbqs_simplified.txt`.

The `traj-compression/` sources used by `benchmark.py` are direct
adaptations of publicly-available code; consult each subdirectory's own
README for the original attribution.

---

## Dataset

The original dataset is the T-Drive taxi trajectory dataset:

> Jing Yuan, Yu Zheng, Xing Xie, and Guangzhong Sun. *Driving with knowledge
> from the physical world.* In KDD '11. ACM.
>
> Jing Yuan, Yu Zheng, Chengyang Zhang, Wenlei Xie, Xing Xie, Guangzhong
> Sun, and Yan Huang. *T-Drive: driving directions based on taxi
> trajectories.* In GIS '10. ACM.

Source on Kaggle: <https://www.kaggle.com/datasets/arashnic/tdriver>

`download_dataset.py` fetches it via `kagglehub`. Only the first 102
trajectories are committed to the repository as `data/<id>/original.txt`;
download the full set yourself if you need more.

---

## Acknowledgements

- The baseline code in `algorithms/` is modified from
  <https://github.com/MingjiHan99/Trajectory-Simplification-Algorithm>.
- The `traj-compression/` baselines (DP, DOTS, SQUISH) are adapted from
  publicly-available implementations of those algorithms.
- The Fréchet distance computation uses
  [FrechetDist.jl](https://github.com/ingomueller-net/FrechetDist.jl).

---

## TODO

- Implementation improvements: reduce the number of `std::vector<Point>` ↔
  `Polygon_2` conversions in `get_points_from_grid`; trim the F/G polygon
  collectors in the streaming loop (currently they grow unbounded and
  cause the multi-GB RSS runaways that the benchmark's memcap exists to
  catch).
- Refactor: collapse the two parallel representations of a point list
  (`data/<id>/original.txt` and `data/taxi/<id>.txt`) into a single
  `N x y` format everywhere.
