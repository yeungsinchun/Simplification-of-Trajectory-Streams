## Simplification of Trajectory Streams (draft)

This repo implements the delta-simplification from the paper [Simplification of Trajectory Streams](https://arxiv.org/abs/2503.23025). It includes a Qt viewer and a CLI to simplify trajectories and optionally compute (approximate) Fréchet distance between original and simplified polylines.

Tested primarily on macOS (arm64). Other platforms may work with the right dependencies.

https://github.com/user-attachments/assets/7593b5d8-ca68-4b7b-af6f-de76cd1c7689

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
  - `taxi/<id>/` — normalized inputs
  - `taxi_simplified/<id>/` — outputs per id: `original.txt`, `simplified.txt`, and `<algo>_simplified.txt`
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
- `main` — baseline algorithms runner (runs all selected algorithms on one id)
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
---

The dataset used is [here](https://www.kaggle.com/datasets/arashnic/tdriver). If you want the full dataset, you will need to download the dataset yourself as only part of the data (the first 102) are included in this repository.

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
- `operba_simplified.txt`
- `fbqs_simplified.txt`

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

- Implementation improvements (e.g., `get_points_from_grid`, reduce the number of `std::vector<Point>` to `Polygon_2` conversion maybe)

- Refactor: duplicated representations of points in taxi/<id>/original.txt and taxi/1/.txt. It's better to use the <N> <x1> <y1> <x2> <y2> representation cuz it's easier to read.