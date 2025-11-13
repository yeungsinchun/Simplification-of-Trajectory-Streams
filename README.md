## Simplification of Trajectory Streams (draft)

This repo implements ideas from the paper [Simplification of Trajectory Streams](https://arxiv.org/abs/2503.23025). It includes a Qt viewer and a CLI to simplify trajectories and optionally compute (approximate) Fréchet distance between original and simplified polylines.

Tested primarily on macOS (arm64). Other platforms may work with the right dependencies.

---

## Prerequisites

- C++ toolchain with C++23 support
- CMake ≥ 3.16
- CGAL
- Qt 6 Widgets (QtCore, QtGui, QtWidgets) for the GUI viewer
- Optional (for frechet distance computation(`--dist`)):
	- Python 3
	- NumPy
    - frechetlib

On macOS with Homebrew you can install them with

```bash
brew install cmake cgal qt@6
pip3 install frechetlib
```

---

## Folder structure

- `simplify.cpp` — main application (CLI + Qt viewer)
- `plot_curves.cpp` — Qt viewer to overlay original + multiple simplified curves with legend and counts
- `frechet` — Python helper to compute (approximate) Fréchet distance from outputs
- `algorithms/` — baseline algorithms we can run for comparison (DP, OPERB, OPERBA, FBQS)
- `tools/normalize.cpp` — C++ tool to normalize raw CSV taxi logs into `data/taxi/<id>.txt`.
- `data/`
  - `taxi/` — normalized inputs as text: `../data/taxi/<id>.txt`
  - `taxi_simplified/<id>/` — per-id outputs: `original.txt`, `simplified.txt`, and `<algo>_simplified.txt`
- `build/` — build directory (create this yourself)

---

## Build and install

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=release ..
cmake --build .
```

Targets that will be built inside `build/`:

- `simplify` — main app and viewer
- `plot_curves` — overlay viewer for all curves for a given id
- `main` — algorithms runner (runs all selected algorithms on one id)
- `normalize` — normalize raw CSV taxi logs into `data/taxi/<id>.txt`

---

## Usage

Run from the `build/` directory after a successful build.

Basic headless run (read input by id and write outputs):

```bash
./simplify --in 1 --out
```

Show GUI viewer while simplifying:

```bash
./simplify --in 1 --gui
```

Full workflow (GUI + write outputs + run distance later when GUI closes):

```bash
./simplify --in 1 --gui --out --dist
```

Tuning parameters at runtime:

```bash
./simplify --in 1 -d 1500 -e 0.5 --gui
```

CLI flags summary (simplify):

- `--in <id>` — read from `../data/taxi/<id>.txt`
- `--out` — write outputs into `../data/taxi_simplified/<id>/`
- `--dist` — compute Fréchet distance by running `../frechet -id <id>` after GUI closes
- `--gui` — show the viewer
- `-d <delta>` — override delta
- `-e <epsilon>` — override epsilon
- `-F`, `-G`, `-S` — show the intermediate structures in the viewer (refer to the paper for description of these structures)
- `-h` — print help

Notes:
- With `--gui`, distance computation is deferred until you close the viewer.

---

### Normalize raw taxi logs

Convert raw CSV logs to normalized XY text files used by this repo:

```bash
cmake --build . --target normalize
./normalize --all          # process all ids found under taxi_log_2008_by_id
# or a single id
./normalize -n 16
```

Behavior:

- Reads the last two comma-separated fields from each CSV line as (x, y).
- Maps X linearly so that min X -> -8000 and max X -> 8000.
- Scales Y linearly around its mean using the same factor as X; if the resulting Y range exceeds [-8000, 8000], it is further linearly scaled (around 0) to fit within [-8000, 8000]. No clipping is performed.
- Writes to `../data/taxi/<id>.txt` with first line N, followed by N lines `x y`.

---

The dataset used is [here](apc01.safelinks.protection.outlook.com/?url=https%3A%2F%2Fwww.kaggle.com%2Fdatasets%2Farashnic%2Ftdriver&data=05%7C02%7Cscyeungaf%40connect.ust.hk%7Ca851043263f44a03421908de04be8c83%7C6c1d415239d044ca88d9b8d6ddca0708%7C1%7C0%7C638953413608899087%7CUnknown%7CTWFpbGZsb3d8eyJFbXB0eU1hcGkiOnRydWUsIlYiOiIwLjAuMDAwMCIsIlAiOiJXaW4zMiIsIkFOIjoiTWFpbCIsIldUIjoyfQ%3D%3D%7C0%7C%7C%7C&sdata=rHbf662%2BB4zP8HyBp6ZFJzSGKxrIMbdAuLFfjdzJmoY%3D&reserved=0). If you want the full dataset, you will need to download the dataset yourself as only part of the data (the first 102) are included in this repository.

---

## Algorithms runner (main)

Run the baseline algorithms on one id and write per-algorithm outputs:

```bash
cmake --build . --target main
./main <id> <error_bound>
# example
./main 1 100
```

This reads `../data/taxi/<id>.txt`, runs all configured algorithms, and writes to `../data/taxi_simplified/<id>/`:

- `dp_simplified.txt`
- `operb_simplified.txt`
- `operba_simplified.txt` (when available)
- `fbqs_simplified.txt`

Each file uses two-line format:

- Line 1: space-separated x values
- Line 2: space-separated y values

Notes:

- The set of algorithms compiled may vary (e.g., OPERBA optional). The tool runs all available ones.

---

## Overlay viewer (plot_curves)

Display original plus selected simplified curves with legend and point counts:

```bash
cmake --build . --target plot_curves
./plot_curves <id> [--all | -dp -operb -operba -fbqs -simplify]

# examples
./plot_curves 1 --all
./plot_curves 1 -dp -operb
./plot_curves 1 -simplify
```

It loads from `../data/taxi_simplified/<id>/*.txt` (and the original curve if available).
Special debug overlays (polygons/points) use distinct styles and appear in the legend with counts.

---

## Fréchet (optional)

The helper `frechet` computes (approximate) distances between the original and available simplified curves for an id:

```bash
./frechet <id>
```

This prints distances for `simplified.txt` and each `<algo>_simplified.txt` found under `../data/taxi_simplified/<id>/`.

---

## Acknowledgements

We are grateful to the following papers for datasets:

1. Jing Yuan, Yu Zheng, Xing Xie, and Guangzhong Sun. Driving with knowledge from the physical world. In KDD ’11. ACM.
2. Jing Yuan, Yu Zheng, Chengyang Zhang, Wenlei Xie, Xing Xie, Guangzhong Sun, and Yan Huang. T-Drive: driving directions based on taxi trajectories. In GIS ’10. ACM.

The code within the `algorithms` folder that we benchmark our algorithm against is modified from
https://github.com/MingjiHan99/Trajectory-Simplification-Algorithm

---

## TODO

- Implementation improvements (e.g., `get_points_from_grid`, reduce the `std::vector<Point>` to `Polygon_2` conversion maybe)