#!/usr/bin/env python3
"""Serial, native core_ms benchmark with measured continuous Frechet bounds.

See README.md "Benchmarking — SOTS vs DOTS / SQUISH / DP" for the timing and calibration contracts.
No process-wall-time fallback, synthetic baseline, or parallel worker is used.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import re
import select
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
ALGORITHMS = ("sots", "dots", "squish", "dp")
EPSILONS = (299.0, 30.0, 5.0, 0.5, 0.1)
BOUNDS = (300.0, 1000.0)
THREAD_ENV = {key: "1" for key in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "JULIA_NUM_THREADS",
    "JULIA_NUM_PRECOMPILE_TASKS")}
MARKERS = {"sots": "SIMPLIFY", "reference": "SIMPLIFY", "dots": "DOTS",
           "dp": "DP", "squish": "SQUISH"}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def curve_info(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0])
    points = [list(map(float, line.split())) for line in lines[1:]]
    if n < 1 or len(points) != n or any(
            len(p) != 2 or not all(map(math.isfinite, p)) for p in points):
        raise ValueError(f"Invalid curve: {path}")
    return n, digest(path)


def core_ms(output, algorithm):
    hits = re.findall(rf"^{MARKERS[algorithm]}_CORE_MS:\s*(\S+)\s*$",
                      output, re.MULTILINE)
    if len(hits) != 1:
        raise ValueError(f"Expected one core timer for {algorithm}: {output}")
    value = float(hits[0])
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"Invalid core timer: {value}")
    return value


def within_bound(distance, bound):
    # Only floating-point slack: e.g. FrechetDist reports 300.00000000000017
    # for a curve exactly on the 300-unit boundary.
    return distance <= bound + max(1e-7, 1e-9 * bound)


class Frechet:
    def __init__(self, env):
        self.process = subprocess.Popen(
            ["julia", "--startup-file=no", str(ROOT / "scripts/fair_frechet.jl")],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, env=env)

    def distance(self, original, simplified):
        if any(c in str(path) for path in (original, simplified) for c in "\t\n"):
            raise ValueError("Curve paths must not contain tabs or newlines")
        self.process.stdin.write(f"{original}\t{simplified}\n")
        self.process.stdin.flush()
        if not select.select([self.process.stdout], [], [], 180)[0]:
            raise TimeoutError("Continuous Frechet worker timed out")
        line = self.process.stdout.readline().strip()
        value = float(line)  # Errors and missing results fail closed.
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"Invalid continuous Frechet distance: {line}")
        return value

    def close(self):
        self.process.stdin.close()
        try:
            self.process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            self.process.kill()
            self.process.wait()
        self.process.stdout.close()


class Runner:
    def __init__(self, build, env, reference=None):
        self.env = env
        self.binaries = {a: build / ("simplify" if a == "sots" else a)
                         for a in ALGORITHMS}
        if reference:
            self.binaries["reference"] = reference
        for path in self.binaries.values():
            if not path.is_file():
                raise FileNotFoundError(f"Build the native CMake target: {path}")

    def run(self, algorithm, dataset, control, epsilon=0.5, delta=200):
        directory = ROOT / "data" / str(dataset)
        original = directory / "original.txt"
        binary = str(self.binaries[algorithm])
        if algorithm in ("sots", "reference"):
            output = directory / "simplify.txt"
            args = [binary, str(dataset), "-e", str(epsilon), "-d", str(delta)]
        elif algorithm == "dots":
            output = directory / "dots_simplified.txt"
            args = [binary, str(dataset), "-lssd", str(control)]
        else:
            output = directory / f"fair_{algorithm}.txt"
            args = [binary, str(original), str(control), str(output)]
        # Never allow an old output to stand in for a failed writer.
        output.unlink(missing_ok=True)
        result = subprocess.run(args, capture_output=True, text=True,
                                env=self.env, cwd=ROOT, timeout=120, check=True)
        ms = core_ms(result.stdout + "\n" + result.stderr, algorithm)
        n, sha = curve_info(output)
        return {"core_ms": ms, "points": n, "sha256": sha}, output


def calibrate(runner, frechet, dataset, bound, algorithm):
    """Verified feasible selection; heuristic search does not assume/prove
    monotone continuous Frechet error or globally minimal retained points."""
    original = ROOT / "data" / str(dataset) / "original.txt"
    n, _ = curve_info(original)
    trials = []
    cache = {}

    def evaluate(control):
        if control not in cache:
            info, path = runner.run(algorithm, dataset, control)
            distance = frechet.distance(original, path)
            row = {"control": control, "points": info["points"],
                   "sha256": info["sha256"], "frechet": distance,
                   "feasible": within_bound(distance, bound)}
            trials.append(row)
            cache[control] = row
        return cache[control]["feasible"]

    if algorithm == "squish":
        # The pinned implementation's buffer<=2 path can omit the last point.
        # Calibrate only its valid buffer>=3 regime (or full retention for n<3).
        low, high = min(3, n), n
        def ratio(k):
            return 1.0 if k == n else (k + 0.25) / n
        evaluate(ratio(high))
        while low < high:
            middle = (low + high) // 2
            if evaluate(ratio(middle)):
                high = middle
            else:
                low = middle + 1
        evaluate(ratio(low))
    else:
        low = 0.0 if algorithm == "dp" else 1e-9
        high = bound if algorithm == "dp" else bound * bound
        if not evaluate(low):
            raise RuntimeError(f"{algorithm} cannot meet bound {bound} on {dataset}")
        for _ in range(20):
            if not evaluate(high):
                break
            low, high = high, high * 2
        for _ in range(14):
            middle = (low + high) / 2
            if evaluate(middle):
                low = middle
            else:
                high = middle
    feasible = [row for row in trials if row["feasible"]]
    if not feasible:
        raise RuntimeError(f"No measured feasible {algorithm} control")
    chosen = min(feasible, key=lambda row: (row["points"], -row["frechet"]))
    return {"dataset": dataset, "bound": bound, "algorithm": algorithm,
            "chosen": chosen, "trials": trials}


def save(path, obj):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def summarize(samples):
    return {"samples_ms": samples, "mean_ms": statistics.mean(samples),
            "std_ms": statistics.stdev(samples), "min_ms": min(samples),
            "max_ms": max(samples)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=("before", "after"), required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--reference", type=Path,
                        help="Preserved pre-optimization simplify binary; required after")
    parser.add_argument("--output", type=Path, default=ROOT / "results/fair-core")
    parser.add_argument("--ids", nargs="+", type=int, default=list(range(11, 31)))
    parser.add_argument("--epsilons", nargs="+", type=float, default=EPSILONS)
    parser.add_argument("--bounds", nargs="+", type=float, default=BOUNDS)
    parser.add_argument("--runs", type=int, default=10)
    args = parser.parse_args()
    if args.runs < 2 or (args.phase == "after" and not args.reference):
        parser.error("At least two runs and an after-phase reference are required")
    if any(x <= 0 or not math.isfinite(x) for x in (*args.epsilons, *args.bounds)):
        parser.error("Epsilons and bounds must be positive and finite")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    destination = output / f"{args.phase}.json"
    if destination.exists():
        parser.error(f"Refusing to overwrite completed evidence: {destination}")
    env = {**os.environ, **THREAD_ENV}
    runner = Runner(args.build_dir.resolve(), env,
                    args.reference.resolve() if args.reference else None)
    started = time.time()
    inputs = {str(i): {"points": curve_info(ROOT / f"data/{i}/original.txt")[0],
                       "sha256": digest(ROOT / f"data/{i}/original.txt")}
              for i in args.ids}
    config = {"ids": args.ids, "epsilons": list(args.epsilons),
              "bounds": list(args.bounds), "runs": args.runs,
              "seed": 20260928, "warmups": 1, "inputs": inputs}
    prior = None
    if args.phase == "after":
        prior = json.loads((output / "before.json").read_text())
        if prior["config"] != config:
            raise ValueError("Before/after input or parameter mismatch")
        for algorithm in ("dots", "squish", "dp"):
            if digest(runner.binaries[algorithm]) != prior["binaries"][algorithm]:
                raise ValueError(f"Changed baseline binary: {algorithm}")
        if digest(runner.binaries["reference"]) != prior["binaries"]["sots"]:
            raise ValueError("Reference is not the measured before binary")
    result = {"phase": args.phase, "config": config, "cases": [],
              "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "platform": platform.platform(), "thread_env": THREAD_ENV,
              "binaries": {a: digest(b) for a, b in runner.binaries.items()},
              "git_head": subprocess.check_output(
                  ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()}
    calibration_path = output / "calibration.json"
    calibrations = json.loads(calibration_path.read_text()) if calibration_path.exists() else {}
    rng = random.Random(config["seed"])
    algorithms = list(runner.binaries)
    frechet = Frechet(env)
    try:
        for dataset in args.ids:
            original = ROOT / f"data/{dataset}/original.txt"
            for bound in args.bounds:
                controls = {}
                for algorithm in ("dots", "squish", "dp"):
                    key = f"{dataset}/{bound:g}/{algorithm}"
                    if key not in calibrations:
                        if prior:
                            raise ValueError(f"Missing frozen calibration: {key}")
                        print(f"Calibrating {key}", flush=True)
                        calibrations[key] = calibrate(runner, frechet, dataset, bound, algorithm)
                        calibrations[key]["input_sha256"] = inputs[str(dataset)]["sha256"]
                        calibrations[key]["binary_sha256"] = digest(runner.binaries[algorithm])
                        save(calibration_path, calibrations)
                    calibration = calibrations[key]
                    if (calibration["input_sha256"] != inputs[str(dataset)]["sha256"] or
                            calibration["binary_sha256"] != digest(runner.binaries[algorithm])):
                        raise ValueError(f"Stale calibration: {key}")
                    controls[algorithm] = calibration["chosen"]["control"]
                for epsilon in args.epsilons:
                    delta = float(f"{bound / (1 + epsilon):.15g}")
                    case = {"dataset": dataset, "epsilon": epsilon, "delta": delta,
                            "bound": bound, "algorithms": {}}
                    samples = {a: [] for a in algorithms}
                    identities = {}
                    # Warmup and timed repetitions use shuffled round-robin order.
                    # Every invocation starts a new process; internal caches are cold.
                    for repetition in range(args.runs + 1):
                        order = algorithms.copy()
                        rng.shuffle(order)
                        for algorithm in order:
                            info, path = runner.run(algorithm, dataset, controls.get(algorithm), epsilon, delta)
                            identity = (info["points"], info["sha256"])
                            if repetition == 0:
                                identities[algorithm] = identity
                                distance = frechet.distance(original, path)
                                if not within_bound(distance, bound):
                                    raise ValueError(f"Bound exceeded: {dataset}/{epsilon}/{bound}/{algorithm}: {distance}")
                                case["algorithms"][algorithm] = {
                                    "control": controls.get(algorithm), "points": info["points"],
                                    "sha256": info["sha256"], "frechet": distance}
                            elif identities[algorithm] != identity:
                                raise ValueError(f"Nondeterministic output: {algorithm}")
                            else:
                                samples[algorithm].append(info["core_ms"])
                    if prior:
                        before_case = next(c for c in prior["cases"] if
                            (c["dataset"], c["epsilon"], c["bound"]) == (dataset, epsilon, bound))
                        for algorithm in ALGORITHMS:
                            now, old = case["algorithms"][algorithm], before_case["algorithms"][algorithm]
                            if any(now[k] != old[k] for k in ("control", "sha256", "frechet", "points")):
                                raise ValueError(f"Output changed before/after: {dataset}/{epsilon}/{bound}/{algorithm}")
                        if identities["reference"] != identities["sots"]:
                            raise ValueError("Reference and candidate curves differ")
                    for algorithm in algorithms:
                        case["algorithms"][algorithm].update(summarize(samples[algorithm]))
                    result["cases"].append(case)
                    save(output / f"{args.phase}.partial.json", result)
                    print(f"{args.phase}: {len(result['cases'])} cases; id={dataset} e={epsilon:g} bound={bound:g}", flush=True)
        result["elapsed_seconds"] = time.time() - started
        save(destination, result)
        (output / f"{args.phase}.partial.json").unlink()
        with (output / f"{args.phase}.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["dataset", "epsilon", "delta", "bound", "algorithm", "control",
                             "points", "frechet", "mean_core_ms", "sample_std_ms", "speedup_vs_sots"])
            for case in result["cases"]:
                for algorithm, row in case["algorithms"].items():
                    speedup = case["algorithms"]["sots"]["mean_ms"] / row["mean_ms"] if row["mean_ms"] else ""
                    writer.writerow([case[k] for k in ("dataset", "epsilon", "delta", "bound")] +
                                    [algorithm, row["control"], row["points"], row["frechet"],
                                     row["mean_ms"], row["std_ms"], speedup])
    finally:
        frechet.close()


if __name__ == "__main__":
    main()
