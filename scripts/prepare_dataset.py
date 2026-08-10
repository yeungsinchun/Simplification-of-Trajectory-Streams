#!/usr/bin/env python3
"""Download and normalize the first N T-Drive trajectories for this project."""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
from pathlib import Path

DATASET = "arashnic/tdriver"
REPO_ROOT = Path(__file__).resolve().parent.parent
RAW_DIRECTORY = REPO_ROOT / "taxi_log_2008_by_id"
OUTPUT_DIRECTORY = REPO_ROOT / "data"


def find_raw_directory(download_root: Path) -> Path | None:
    """Find the expected raw trajectory directory in a Kaggle download."""
    candidates = [
        download_root / "taxi_log_2008_by_id",
        download_root / "release" / "taxi_log_2008_by_id",
    ]
    candidates.extend(download_root.rglob("taxi_log_2008_by_id"))
    return next((candidate for candidate in candidates if candidate.is_dir()), None)


def ensure_raw_data() -> Path:
    """Use local raw data when present; otherwise download it through Kagglehub."""
    if RAW_DIRECTORY.is_dir():
        return RAW_DIRECTORY

    try:
        import kagglehub
    except ImportError as error:
        raise RuntimeError(
            "Missing dependency: install it with `python3 -m pip install kagglehub`."
        ) from error

    try:
        download_root = Path(kagglehub.dataset_download(DATASET))
    except Exception as error:
        raise RuntimeError(
            "Unable to download T-Drive. Configure Kaggle credentials and retry.\n"
            f"Kagglehub error: {error}"
        ) from error

    source = find_raw_directory(download_root)
    if source is None:
        raise RuntimeError(
            f"Downloaded dataset did not contain taxi_log_2008_by_id under {download_root}."
        )

    RAW_DIRECTORY.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(source, RAW_DIRECTORY)
    return RAW_DIRECTORY


def parse_last_two_csv_fields(line: str) -> tuple[float, float] | None:
    """Read x and y from the final two comma-separated fields."""
    try:
        fields = next(csv.reader([line]))
        if len(fields) < 2:
            return None
        return float(fields[-2].strip()), float(fields[-1].strip())
    except (ValueError, csv.Error):
        return None


def normalize(raw_path: Path) -> list[tuple[float, float]]:
    points = []
    for line in raw_path.read_text().splitlines():
        if line.strip() and (point := parse_last_two_csv_fields(line)) is not None:
            points.append(point)
    if not points:
        raise RuntimeError(f"No valid points parsed from {raw_path}")

    xs, ys = zip(*points)
    min_x = min(xs)
    x_range = max(xs) - min_x
    y_mean = sum(ys) / len(ys)
    max_y_deviation = max(abs(y - y_mean) for y in ys)
    x_scale = 20000.0 / x_range if x_range > 0 else 1.0
    y_scale = 10000.0 / max_y_deviation if max_y_deviation > 0 else 1.0
    scale = min(x_scale, y_scale)
    return [(-10000.0 + (x - min_x) * scale, (y - y_mean) * scale) for x, y in points]


def write_curve(trajectory_id: int, points: list[tuple[float, float]]) -> None:
    output_path = OUTPUT_DIRECTORY / str(trajectory_id) / "original.txt"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as output:
        output.write(f"{len(points)}\n")
        for x, y in points:
            output.write(f"{x:.15g} {y:.15g}\n")
    print(f"Wrote {output_path} ({len(points)} points)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download T-Drive when needed and normalize trajectories 1 through N."
    )
    parser.add_argument(
        "--size", type=int, required=True,
        help="number of trajectories to prepare, starting at ID 1",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.size < 1:
        print("--size must be at least 1.", file=sys.stderr)
        return 2

    try:
        raw_directory = ensure_raw_data()
        missing = [str(trajectory_id) for trajectory_id in range(1, args.size + 1)
                   if not (raw_directory / f"{trajectory_id}.txt").is_file()]
        if missing:
            raise RuntimeError(
                f"Raw dataset is missing required trajectory IDs: {', '.join(missing)}"
            )
        for trajectory_id in range(1, args.size + 1):
            write_curve(trajectory_id, normalize(raw_directory / f"{trajectory_id}.txt"))
    except RuntimeError as error:
        print(error, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
