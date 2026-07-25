#!/usr/bin/env python3
"""Benchmark simplify against DOTS using algorithm-only wall times."""

from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import statistics
import subprocess
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SIMPLIFY = REPO / "release" / "simplify"
DOTS = REPO / "release" / "DOTS"
FRECHET = REPO / "scripts" / "frechet"
DATA = REPO / "data"
OUTPUT = REPO / "results" / "compare_dots_core.csv"
DOTS_THRESHOLD = 1000.0
EPSILONS = (0.5, 0.666, 0.8, 1.0)
DELTA_SCALES = (0.75, 1.0, 1.25, 1.5)
REFINE_SCALES = (0.9, 0.95, 1.05, 1.1)
MATCH_TOLERANCE = 0.05
CORE_REPEATS = 20

DOTS_TIME = re.compile(r"DOTS_CORE_MS\s+([0-9]+(?:\.[0-9]+)?)")
FRECHET_LINE = re.compile(r"[^:]+:\s*([-+0-9.eE]+)\s*$")


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, cwd=REPO, capture_output=True, text=True, check=False)


def point_count(path: Path) -> int:
    with path.open() as stream:
        value = int(stream.readline().strip())
    if value <= 0:
        raise RuntimeError(f"empty output: {path}")
    return value


def frechet_batch(idv: int, paths: list[Path]) -> dict[Path, float]:
    command = [str(FRECHET), "--id", str(idv)]
    unique_paths = list(dict.fromkeys(paths))
    command.extend(arg for path in unique_paths for arg in ("--batch", str(path)))
    result = run(command)
    if result.returncode != 0:
        raise RuntimeError(f"Frechet failed for ID {idv}: {result.stderr.strip()}")
    values: dict[Path, float] = {}
    by_name = {path.name: path for path in unique_paths}
    for line in result.stdout.splitlines():
        match = FRECHET_LINE.fullmatch(line.strip())
        if match:
            name = line.split(":", 1)[0].strip()
            if name in by_name:
                values[by_name[name]] = float(match.group(1))
    if len(values) != len(unique_paths):
        raise RuntimeError(f"Frechet returned {len(values)}/{len(unique_paths)} results for ID {idv}")
    return values


def run_dots(idv: int, repeats: int, output: Path) -> tuple[list[float], int, Path]:
    command = [str(DOTS), str(idv), str(DOTS_THRESHOLD)]
    times: list[float] = []
    for _ in range(repeats + 1):
        result = run(command)
        if result.returncode != 0:
            raise RuntimeError(f"DOTS failed for ID {idv}: {result.stderr.strip()}")
        matches = DOTS_TIME.findall(result.stderr)
        if len(matches) != 1:
            raise RuntimeError(f"DOTS_CORE_MS missing or duplicated for ID {idv}")
        times.append(float(matches[0]))
    # The first run is a warm-up; all later runs are measured.
    times = times[1:]
    source = DATA / str(idv) / "dots_simplified.txt"
    if not source.exists():
        raise RuntimeError(f"DOTS output missing: {source}")
    shutil.copy2(source, output)
    return times, point_count(output), output


def run_simplify_once(idv: int, epsilon: float, delta: float, output: Path) -> float:
    start = time.monotonic()
    result = run([str(SIMPLIFY), str(idv), "--out", "-e", str(epsilon), "-d", str(delta)])
    elapsed_ms = (time.monotonic() - start) * 1000.0
    if result.returncode != 0:
        raise RuntimeError(f"simplify failed for ID {idv}, e={epsilon}, d={delta}: {result.stderr.strip()}")
    source = DATA / str(idv) / "simplify.txt"
    if not source.exists():
        raise RuntimeError(f"simplify output missing: {source}")
    shutil.copy2(source, output)
    return elapsed_ms


def make_candidates(idv: int, work: Path, parameters: list[tuple[float, float]]) -> list[tuple[float, float, Path, int]]:
    candidates: list[tuple[float, float, Path, int]] = []
    for epsilon, delta in parameters:
        path = work / f"simplify_e{epsilon}_d{delta:.8g}.txt"
        try:
            run_simplify_once(idv, epsilon, delta, path)
            candidates.append((epsilon, delta, path, point_count(path)))
        except (RuntimeError, ValueError) as exc:
            print(f"  calibration failed e={epsilon} d={delta}: {exc}")
    return candidates


def calibrate_simplify(idv: int, work: Path, dots_path: Path) -> tuple[float, float, Path, int, float, float]:
    dots_error = frechet_batch(idv, [dots_path])[dots_path]
    if dots_error == 0.0:
        parameters = [(epsilon, 0.0) for epsilon in EPSILONS]
    else:
        parameters = [
            (epsilon, dots_error / (1.0 + epsilon) * scale)
            for epsilon in EPSILONS
            for scale in DELTA_SCALES
        ]
    candidates = make_candidates(idv, work, parameters)
    if not candidates:
        raise RuntimeError(f"no simplify calibration candidate succeeded for ID {idv}")
    errors = frechet_batch(idv, [candidate[2] for candidate in candidates])
    scored = [(*candidate, errors[candidate[2]]) for candidate in candidates]
    chosen = min(scored, key=lambda item: abs(item[4] - dots_error) / max(dots_error, 1e-12))
    mismatch = abs(chosen[4] - dots_error) / max(dots_error, 1e-12)
    if mismatch > MATCH_TOLERANCE:
        epsilon, delta = chosen[0], chosen[1]
        refinements = make_candidates(idv, work, [(epsilon, delta * scale) for scale in REFINE_SCALES])
        if refinements:
            refine_errors = frechet_batch(idv, [candidate[2] for candidate in refinements])
            scored.extend((*candidate, refine_errors[candidate[2]]) for candidate in refinements)
            chosen = min(scored, key=lambda item: abs(item[4] - dots_error) / max(dots_error, 1e-12))
    return (*chosen, dots_error)


def measure_simplify(idv: int, epsilon: float, delta: float, repeats: int, output: Path) -> list[float]:
    times: list[float] = []
    for _ in range(repeats + 1):
        times.append(run_simplify_once(idv, epsilon, delta, output))
    return times[1:]


def iqr(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    return statistics.quantiles(values, n=4, method="inclusive")[2] - statistics.quantiles(values, n=4, method="inclusive")[0]


def finite(value: float) -> bool:
    return math.isfinite(value) and value >= 0.0


def write_rows(output: Path, rows: list[dict[str, object]]) -> None:
    fields = list(rows[0]) if rows else []
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def benchmark(ids: range, repeats: int, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for idv in ids:
        print(f"ID {idv}: DOTS threshold={DOTS_THRESHOLD:g}", flush=True)
        work = DATA / str(idv) / "core_benchmark"
        work.mkdir(parents=True, exist_ok=True)
        dots_path = work / "DOTS.txt"
        dots_times, dots_points, dots_result_path = run_dots(idv, repeats, dots_path)
        epsilon, delta, _calibration_path, simplify_points, simplify_error, dots_error = calibrate_simplify(
            idv, work, dots_result_path
        )
        simplify_path = work / "simplify.txt"
        simplify_times = measure_simplify(idv, epsilon, delta, repeats, simplify_path)
        rel_mismatch = abs(simplify_error - dots_error) / max(dots_error, 1e-12)
        quality = "ok" if rel_mismatch <= MATCH_TOLERANCE else "mismatch"
        row = {
            "id": idv,
            "dots_threshold": DOTS_THRESHOLD,
            "dots_points": dots_points,
            "simplify_points": simplify_points,
            "dots_frechet": dots_error,
            "simplify_frechet": simplify_error,
            "relative_frechet_mismatch": rel_mismatch,
            "quality_status": quality,
            "simplify_epsilon": epsilon,
            "simplify_delta": delta,
            "dots_median_ms": statistics.median(dots_times),
            "dots_iqr_ms": iqr(dots_times),
            "simplify_median_ms": statistics.median(simplify_times),
            "simplify_iqr_ms": iqr(simplify_times),
            "dots_runs": repeats,
            "simplify_runs": repeats,
        }
        if not all(finite(float(row[key])) for key in ("dots_median_ms", "dots_iqr_ms", "simplify_median_ms", "simplify_iqr_ms")):
            raise RuntimeError(f"non-finite timing for ID {idv}")
        rows.append(row)
        write_rows(output, rows)
        print(
            f"  DOTS {row['dots_median_ms']:.3f} ms, simplify {row['simplify_median_ms']:.3f} ms, "
            f"Fréchet {dots_error:.3f}/{simplify_error:.3f}, {quality}", flush=True
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=int, default=1)
    parser.add_argument("--end", type=int, default=10)
    parser.add_argument("--repeats", type=int, default=CORE_REPEATS)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if args.start > args.end or args.repeats < 1:
        parser.error("require start <= end and repeats >= 1")
    if not SIMPLIFY.is_file() or not DOTS.is_file() or not FRECHET.is_file():
        parser.error(f"missing benchmark binary or Frechet executable under {REPO}")
    benchmark(range(args.start, args.end + 1), args.repeats, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
