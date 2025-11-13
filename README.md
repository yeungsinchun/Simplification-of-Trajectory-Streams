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
- `data/`
	- `taxi/` — raw inputs as text: `../data/taxi/<id>.txt` (see format below)
	- `taxi_simplified/<id>/` — outputs written by `--out`: `original.txt`, `simplified.txt`
- `script/`
	- `clean` — shell script to pre-process taxi logs into the expected input format
- `frechet` — helper to compute (approximate) Fréchet distance from outputs
- `build/` — out-of-source build directory (create this yourself)
- `README.md` — this file

---

## Build and install

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```

This produces an executable named `simplify` in `build/`.

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

CLI flags summary:

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

## Cleaning and preparing data

There is a helper script to format raw taxi logs to the expected input layout:

```bash
./script/clean --all
```

This enumerates input files numerically by basename (e.g., `1.txt`, `2.txt`, …), processes each one, and writes cleaned results. Check the script header/comments for more options. The dataset used is [here](apc01.safelinks.protection.outlook.com/?url=https%3A%2F%2Fwww.kaggle.com%2Fdatasets%2Farashnic%2Ftdriver&data=05%7C02%7Cscyeungaf%40connect.ust.hk%7Ca851043263f44a03421908de04be8c83%7C6c1d415239d044ca88d9b8d6ddca0708%7C1%7C0%7C638953413608899087%7CUnknown%7CTWFpbGZsb3d8eyJFbXB0eU1hcGkiOnRydWUsIlYiOiIwLjAuMDAwMCIsIlAiOiJXaW4zMiIsIkFOIjoiTWFpbCIsIldUIjoyfQ%3D%3D%7C0%7C%7C%7C&sdata=rHbf662%2BB4zP8HyBp6ZFJzSGKxrIMbdAuLFfjdzJmoY%3D&reserved=0). If you want the full dataset, you will need to download the dataset yourself as only part of the data (the first 102) are included in this repository.

---

## Acknowledgements

We are grateful to the following papers for datasets:

1. Jing Yuan, Yu Zheng, Xing Xie, and Guangzhong Sun. Driving with knowledge from the physical world. In KDD ’11. ACM.
2. Jing Yuan, Yu Zheng, Chengyang Zhang, Wenlei Xie, Xing Xie, Guangzhong Sun, and Yan Huang. T-Drive: driving directions based on taxi trajectories. In GIS ’10. ACM.

The code within the `algorithms` folder that we benchmark our algorithm against is modified from
https://github.com/MingjiHan99/Trajectory-Simplification-Algorithm

---

## TODO

- Implementation improvements (e.g., `get_points_from_grid`)