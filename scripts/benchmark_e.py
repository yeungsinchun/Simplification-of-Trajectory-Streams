#!/usr/bin/env python3
"""
Sweep epsilon to analyse its effect on simplification quality and compression.

For each trajectory ID and each epsilon value:
  - delta = 300 / (1 + epsilon)   (keeps the corridor area roughly constant)
  - Run simplify --in <id> -e <epsilon> -d <delta>
  - Record original_points, simplified_points, epsilon, delta, runtime (SIMPLIFY_CORE_MS)

Results are written to results/benchmark_e.csv.

Usage:
  python scripts/benchmark_e.py [--a 1] [--b 1000] [--workers 8] [--resume]
"""

import argparse
import concurrent.futures
import csv
import functools
import multiprocessing as _mp
import os
import re
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Optional, Tuple

import psutil

REPO_ROOT = Path(__file__).resolve().parent.parent
SIMPLIFY   = REPO_ROOT / "build-release" / "simplify"
DATA_DIR   = REPO_ROOT / "data"
OUTPUT_CSV = REPO_ROOT / "results" / "benchmark_e.csv"

EPSILONS = [0.25, 0.5, 1, 2, 5, 10, 50, 100]

DEFAULT_MEM_CAP_MB  = 1500
DEFAULT_TIMEOUT_S   = 120
MEM_POLL_INTERVAL_S = 0.5

CSV_COLS = [
    "id", "epsilon", "delta",
    "original_points", "simplified_points",
    "core_ms", "status",
]


# ---------------------------------------------------------------------------
# Subprocess helpers (trimmed from benchmark.py)
# ---------------------------------------------------------------------------

def _rss_mb(pid: int) -> Optional[float]:
    try:
        return psutil.Process(pid).memory_info().rss / (1024 * 1024)
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        return None


def _terminate_proc(proc: subprocess.Popen) -> None:
    if proc.poll() is not None:
        return
    try:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            proc.kill()
    except Exception:
        pass
    try:
        proc.wait(timeout=2)
    except subprocess.TimeoutExpired:
        pass


def _watch_mem(proc, mem_cap_mb, tag, state, done_event):
    if mem_cap_mb <= 0:
        return
    cap = mem_cap_mb * 1024 * 1024
    while True:
        if proc.poll() is not None or done_event.is_set():
            return
        rss = _rss_mb(proc.pid)
        if rss is not None and rss * 1024 * 1024 >= cap:
            state["status"] = "memcap"
            print(f"[{tag}] MEMCAP {rss:.0f} MB >= {mem_cap_mb} MB; killing",
                  file=sys.stderr, flush=True)
            _terminate_proc(proc)
            return
        if done_event.wait(MEM_POLL_INTERVAL_S):
            return


def run_cmd(cmd, cwd=None, timeout=120, mem_cap_mb=0, tag="cmd"):
    """Return (status, stdout, stderr). status: 'ok'|'timeout'|'memcap'|'error'."""
    state: dict = {"status": "error"}
    try:
        proc = subprocess.Popen(
            cmd, cwd=cwd, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            start_new_session=True,
        )
    except FileNotFoundError as e:
        return "error", "", f"executable not found: {e}"

    done_event = threading.Event()
    watcher = None
    if mem_cap_mb > 0:
        watcher = threading.Thread(
            target=_watch_mem,
            args=(proc, mem_cap_mb, tag, state, done_event),
            daemon=True,
        )
        watcher.start()

    timed_out = False
    try:
        stdout, stderr = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        print(f"[{tag}] TIMEOUT after {timeout}s", file=sys.stderr, flush=True)
        _terminate_proc(proc)
        try:
            stdout, stderr = proc.communicate(timeout=2)
        except Exception:
            stdout, stderr = "", ""
    finally:
        done_event.set()
        if watcher is not None:
            watcher.join(timeout=2)

    if state.get("status") == "memcap":
        return "memcap", stdout or "", stderr or ""
    if timed_out:
        return "timeout", stdout or "", stderr or ""
    if proc.returncode == 0:
        return "ok", stdout, stderr
    return "error", stdout, stderr


# ---------------------------------------------------------------------------
# Per-(id, epsilon) logic
# ---------------------------------------------------------------------------

def count_points(path: Path) -> Optional[int]:
    try:
        with path.open("r") as f:
            header = f.readline()
            return int(header.strip()) if header else None
    except Exception:
        return None


_CORE_MS_RE = re.compile(r"SIMPLIFY_CORE_MS:\s*([\d.]+)")


def run_one(idv: int, epsilon: float, delta: float,
            timeout: int, mem_cap_mb: int) -> dict:
    """Run simplify for one (id, epsilon) and return a result dict."""
    tag = f"simplify id={idv} e={epsilon} d={delta:.4f}"
    cmd = [
        str(SIMPLIFY), str(idv),
        "-e", str(epsilon),
        "-d", str(delta),
    ]
    t0 = time.monotonic()
    status, stdout, _stderr = run_cmd(cmd, cwd=REPO_ROOT,
                                      timeout=timeout,
                                      mem_cap_mb=mem_cap_mb,
                                      tag=tag)
    elapsed_ms = (time.monotonic() - t0) * 1000

    core_ms: Optional[float] = None
    simplified_points: Optional[int] = None

    if status == "ok":
        m = _CORE_MS_RE.search(stdout)
        if m:
            core_ms = float(m.group(1))
        # simplify writes output to data/<id>/simplify.txt when --out is passed;
        # without --out it still prints "Simplified to N points (X%)." on stdout.
        m2 = re.search(r"Simplified to (\d+) points", stdout)
        if m2:
            simplified_points = int(m2.group(1))

    original_path = DATA_DIR / str(idv) / "original.txt"
    original_points = count_points(original_path)

    return {
        "id": idv,
        "epsilon": epsilon,
        "delta": round(delta, 6),
        "original_points": original_points if original_points is not None else -1,
        "simplified_points": simplified_points if simplified_points is not None else -1,
        "core_ms": round(core_ms, 4) if core_ms is not None else -1,
        "status": status,
    }


def process_id(idv: int, epsilons, timeout, mem_cap_mb,
               skip_keys: Optional[set] = None) -> list:
    original = DATA_DIR / str(idv) / "original.txt"
    if not original.exists():
        return []

    rows = []
    for eps in epsilons:
        key = (idv, eps)
        if skip_keys and key in skip_keys:
            continue
        delta = 300.0 / (1.0 + eps)
        result = run_one(idv, eps, delta, timeout=timeout, mem_cap_mb=mem_cap_mb)
        rows.append(result)
        print(
            f"[ID:{idv}] e={eps:6.2f} d={delta:8.4f} → "
            f"{result['simplified_points']:>6} pts  {result['core_ms']:>10.2f} ms  [{result['status']}]"
        )
    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _process_id_worker(idv, timeout, mem_cap_mb, skip_keys):
    """Top-level worker so ProcessPoolExecutor can pickle it (closures can't)."""
    return process_id(idv, EPSILONS, timeout=timeout,
                      mem_cap_mb=mem_cap_mb, skip_keys=skip_keys)


def _load_done_keys(csv_path: Path) -> set:
    """Return {(id, epsilon)} pairs that already have a valid row (status=='ok')."""
    done = set()
    if not csv_path.exists():
        return done
    try:
        with csv_path.open("r", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    if row.get("status") == "ok":
                        done.add((int(row["id"]), float(row["epsilon"])))
                except (KeyError, ValueError):
                    pass
    except Exception:
        pass
    return done


def main():
    parser = argparse.ArgumentParser(description="Sweep epsilon effect on simplification")
    parser.add_argument("--a", type=int, default=1, help="start id (inclusive, default 1)")
    parser.add_argument("--b", type=int, default=1000, help="end id (inclusive, default 1000)")
    parser.add_argument("--workers", "-w", type=int, default=8,
                        help="parallel worker processes (default 8)")
    parser.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT_S,
                        help=f"per-invocation timeout in seconds (default {DEFAULT_TIMEOUT_S})")
    parser.add_argument("--mem-cap-mb", type=int, default=DEFAULT_MEM_CAP_MB,
                        help=f"RSS cap per child in MB (default {DEFAULT_MEM_CAP_MB}; 0 disables)")
    parser.add_argument("--resume", action="store_true",
                        help="skip (id, epsilon) pairs that already have status=='ok' in the CSV")
    args = parser.parse_args()

    ids = list(range(args.a, args.b + 1))

    skip_keys: set = set()
    if args.resume and OUTPUT_CSV.exists():
        skip_keys = _load_done_keys(OUTPUT_CSV)
        print(f"Resuming: {len(skip_keys)} (id, epsilon) pairs already done, skipping them")

    write_header = not OUTPUT_CSV.exists()
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    _worker = functools.partial(
        _process_id_worker,
        timeout=args.timeout,
        mem_cap_mb=args.mem_cap_mb,
        skip_keys=skip_keys,
    )

    total_rows = 0
    completed = 0
    n = len(ids)

    def _write(rows):
        nonlocal total_rows, write_header
        if not rows:
            return
        with OUTPUT_CSV.open("a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=CSV_COLS)
            if write_header:
                writer.writeheader()
                write_header = False
            writer.writerows(rows)
        total_rows += len(rows)

    print(f"Running epsilon sweep over IDs {args.a}–{args.b} ({n} IDs), "
          f"{len(EPSILONS)} epsilons each, {args.workers} workers")
    print(f"Epsilons: {EPSILONS}")
    print(f"delta = 300 / (1 + epsilon)")

    if args.workers == 1:
        for idv in ids:
            rows = _worker(idv)
            _write(rows)
            completed += 1
            if completed % 50 == 0 or completed == n:
                print(f"Progress: {completed}/{n} IDs ({total_rows} rows)")
    else:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=args.workers,
            mp_context=_mp.get_context("fork"),
        ) as exe:
            future_to_id = {exe.submit(_worker, idv): idv for idv in ids}
            for fut in concurrent.futures.as_completed(future_to_id):
                idv = future_to_id[fut]
                try:
                    rows = fut.result()
                except Exception as e:
                    print(f"[ID:{idv}] worker error: {e}", file=sys.stderr)
                    rows = []
                _write(rows)
                completed += 1
                if completed % 50 == 0 or completed == n:
                    print(f"Progress: {completed}/{n} IDs ({total_rows} rows)")

    print(f"Done. {total_rows} rows written to {OUTPUT_CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
