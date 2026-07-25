#!/usr/bin/env python3
"""
Benchmark script to compare simplify against DP, DOTS, and SQUISH.

For each baseline algorithm:
1. Run the algorithm on original.txt
2. Compute Frechet distance (original vs algo output)
3. For e in {0.5, 0.666, 1.0}:
   - d = frechet_dist / (1 + e)
   - Run simplify with delta=d, epsilon=e
   - Save output to data/<id>/<BASELINE>_against_simplify.txt
   - Compute Frechet distance for simplify output
4. Append results to compare_points.csv

The DP, DOTS and SQUISH binaries are built into release/ by the top-level
CMakeLists.txt (alongside `simplify`), so this script only needs the
release/ directory to exist.
"""
import argparse
import concurrent.futures
import csv
import functools
import os
import re
import statistics
import subprocess
import sys
import threading
import time
import shutil
import signal
from pathlib import Path
from typing import Optional, Tuple

import psutil

REPO_ROOT = Path(__file__).resolve().parent.parent
DP_TOOL = REPO_ROOT / "release" / "DP"
SIMPLIFY = REPO_ROOT / "./release" / "simplify"
DOTS_BIN = REPO_ROOT / "release" / "DOTS"
SQUISH_BIN = REPO_ROOT / "release" / "SQUISH"
FRECHET = REPO_ROOT / "scripts" / "frechet"
DATA_DIR = REPO_ROOT / "data"
EPS_VALUES = [0.5, 0.666, 1.0]
SQUISH_BUFFER = 0.1
OUTPUT_CSV = REPO_ROOT / "compare_points.csv"

# Per-invocation resource limits. These defaults are sized to keep the
# benchmark overnight-stable on a typical workstation: simplify has been
# observed to balloon to multi-GB on pathological (id, delta, epsilon)
# combinations because the algorithm's global F/G polygon collectors and
# CGAL exact-arithmetic Polygon types are never trimmed across iterations.
# A 1.5 GB cap kills such runaways before they OOM the box; the failure is
# recorded as -1 in the CSV and the loop continues with the next ID.
DEFAULT_MEM_CAP_MB = 1500
DEFAULT_ALGO_TIMEOUT_S = 30
DEFAULT_SH_TIMEOUT_S = 300
DEFAULT_FRECHET_TIMEOUT_S = 60
MEM_POLL_INTERVAL_S = 0.5

# Use the juliaup-managed binary directly to avoid the wrapper needing to
# create a lockfile at ~/.julia/juliaup/.juliaup-lock (which is blocked when
# the benchmark runs in a sandbox that forbids writes outside the workspace).
# Falls back to plain "julia" if the versioned binary is missing.
_JULIA_VERSIONED = (
    Path.home() / ".julia" / "juliaup" / "julia-1.12.6+0.aarch64.apple.darwin14"
    / "Julia-1.12.app" / "Contents" / "Resources" / "julia" / "bin" / "julia"
)
JULIA_BIN = str(_JULIA_VERSIONED) if _JULIA_VERSIONED.exists() else "julia"


def _terminate_proc(proc: subprocess.Popen) -> None:
    """Terminate a child process and its descendants. SIGKILL after a grace period."""
    if proc.poll() is not None:
        return
    try:
        # Kill the whole process group so any children (CGAL sub-work, etc.) die too.
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            proc.kill()
    except Exception:
        pass
    try:
        proc.wait(timeout=2)
    except subprocess.TimeoutExpired:
        try:
            proc.kill()
            proc.wait(timeout=1)
        except Exception:
            pass


def _rss_mb(pid: int) -> Optional[float]:
    """Total RSS in MB for the process itself, or None on error.

    We intentionally do NOT walk children(recursive=True): on macOS that
    requires enumerating the global process table (sysctl KERN_PROC_ALL),
    which is blocked in some sandboxes. The Popen'd child is what we care
    about — algorithm binaries don't typically fork worker subprocesses.
    """
    try:
        p = psutil.Process(pid)
        return p.memory_info().rss / (1024.0 * 1024.0)
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        return None


def _watch_mem(proc: subprocess.Popen, mem_cap_mb: int, tag: str,
               state: dict, done_event: threading.Event) -> None:
    """Poll the child's RSS and kill it if it exceeds mem_cap_mb.

    Sets state['status'] to 'memcap' on kill. Runs until either the child
    exits or `done_event` fires. The event-based wakeup avoids the
    `time.sleep(MEM_POLL_INTERVAL_S)` floor in run_cmd's `watcher.join()`:
    if the child has already finished, the parent thread doesn't block
    on the join — it wakes up on the next wait() and exits.
    """
    if mem_cap_mb <= 0:
        return
    cap_bytes = mem_cap_mb * 1024 * 1024
    while True:
        if proc.poll() is not None:
            return
        if done_event.is_set():
            return
        try:
            rss = _rss_mb(proc.pid)
            if rss is not None and rss * 1024 * 1024 >= cap_bytes:
                state['status'] = 'memcap'
                state['peak_mb'] = rss
                print(f"[{tag}] MEMCAP exceeded ({rss:.0f} MB >= {mem_cap_mb} MB); killing PID {proc.pid}",
                      file=sys.stderr, flush=True)
                _terminate_proc(proc)
                return
        except Exception as e:
            # Don't let a flaky psutil call kill the watcher.
            print(f"[{tag}] watcher error: {e}", file=sys.stderr, flush=True)
        # Sleep with an event-based wait so we can be released early when
        # run_cmd() signals the child has exited. A bare time.sleep() here
        # would force run_cmd to block on watcher.join() for up to one
        # full poll interval after the child died — that was inflating
        # every recorded `baseline_time` to ~0.5s on fast algorithms.
        if done_event.wait(MEM_POLL_INTERVAL_S):
            return


def run_cmd(cmd, cwd=None, timeout=120, mem_cap_mb=0, tag="cmd"):
    """Run command and return (status, stdout, stderr).

    status is one of: 'ok', 'timeout', 'memcap', 'error'.
    A background watcher enforces mem_cap_mb (RSS tree, in MB); 0 disables it.
    """
    state: dict = {'status': 'error', 'peak_mb': 0.0}
    try:
        # New session so we can kill the whole process group on overrun.
        proc = subprocess.Popen(
            cmd, cwd=cwd, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            start_new_session=True,
        )
    except FileNotFoundError as e:
        return 'error', '', f"executable not found: {e}"

    watcher = None
    done_event = threading.Event()
    if mem_cap_mb > 0:
        watcher = threading.Thread(
            target=_watch_mem, args=(proc, mem_cap_mb, tag, state, done_event),
            daemon=True
        )
        watcher.start()

    timed_out = False
    try:
        stdout, stderr = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        print(f"[{tag}] TIMEOUT after {timeout}s; killing PID {proc.pid}",
              file=sys.stderr, flush=True)
        _terminate_proc(proc)
        try:
            stdout, stderr = proc.communicate(timeout=2)
        except Exception:
            stdout, stderr = '', ''
    finally:
        # Signal the watcher so it doesn't block on its next sleep cycle.
        done_event.set()
        if watcher is not None:
            watcher.join(timeout=2)

    if state.get('status') == 'memcap':
        return 'memcap', stdout or '', stderr or ''
    if timed_out:
        return 'timeout', stdout or '', stderr or ''
    if proc.returncode == 0:
        return 'ok', stdout, stderr
    return 'error', stdout, stderr


def count_points(path: Path) -> Optional[int]:
    try:
        with path.open('r') as f:
            header = f.readline()
            if not header:
                return None
            return int(header.strip())
    except Exception:
        return None


def _extract_id(path: Path) -> Optional[int]:
    """Extract the trajectory ID from a data path like data/1/original.txt."""
    for p in path.parts:
        if p.isdigit():
            return int(p)
    return None


def frechet(orig_path: Path, simp_path: Path, timeout=120, mem_cap_mb=0,
            tag="frechet") -> Tuple[Optional[float], float, str]:
    """Compute Frechet distance via the frechet Julia script.

    Returns (distance_or_None, elapsed_seconds, status) where status is one of
    'ok', 'timeout', 'memcap', 'error'.
    """
    t0 = time.monotonic()
    traj_id = _extract_id(orig_path)
    if traj_id is None:
        return None, time.monotonic() - t0, 'error'

    cmd = [JULIA_BIN, str(FRECHET), "--in", str(traj_id),
           "--path", str(simp_path), "--raw"]
    status, stdout, _stderr = run_cmd(cmd, cwd=REPO_ROOT, timeout=timeout,
                                     mem_cap_mb=mem_cap_mb, tag=tag)
    elapsed = time.monotonic() - t0
    if status != 'ok':
        return None, elapsed, status
    try:
        return float(stdout.strip()), elapsed, status
    except ValueError:
        return None, elapsed, 'error'


def frechet_batch(idv: int, paths: list,
                  timeout=120, mem_cap_mb=0,
                  tag="frechet_batch") -> Tuple[list, float, str]:
    """Compute Frechet distances for multiple paths in a single Julia invocation.

    Passes all paths as repeated --batch flags so the Julia program can
    compute Frechet against them in a single session. Returns
    ([(path, distance_or_None, status), ...], total_elapsed_seconds,
    overall_status) where overall_status is 'ok' if >=1 succeeded, 'error' if
    all failed, or 'timeout'/'memcap' if the Julia process itself was killed.
    """
    t0 = time.monotonic()

    cmd = [JULIA_BIN, str(FRECHET), "--id", str(idv)]
    for p in paths:
        cmd += ["--batch", str(p)]
    status, stdout, stderr = run_cmd(cmd, cwd=REPO_ROOT, timeout=timeout,
                                     mem_cap_mb=mem_cap_mb, tag=tag)

    elapsed = time.monotonic() - t0

    results = []
    ok_count = 0

    if status != 'ok':
        # Entire Julia call failed (timeout / OOM / crash)
        for p in paths:
            results.append((str(p), None, status))
        return results, elapsed, status

    for line in stdout.strip().splitlines():
        line = line.strip()
        if not line or line.startswith("MISSING") or line.startswith("ERROR"):
            continue
        # Output format: "basename: <number>" (always, even in --raw mode)
        parts = line.rsplit(None, 1)
        if len(parts) == 2 and parts[1][0].isdigit():
            path_hint, dist_str = parts
            # Match by basename suffix
            for p in paths:
                if str(p).endswith(path_hint):
                    try:
                        results.append((str(p), float(dist_str), 'ok'))
                        ok_count += 1
                    except ValueError:
                        results.append((str(p), None, 'error'))
                    break
            else:
                # No basename match; fall back to index order
                idx = len(results)
                if idx < len(paths):
                    try:
                        results.append((str(paths[idx]), float(dist_str), 'ok'))
                        ok_count += 1
                    except ValueError:
                        results.append((str(paths[idx]), None, 'error'))

    overall = 'ok' if ok_count > 0 else 'error'
    return results, elapsed, overall


def _fmt_eps(e: float) -> str:
    """Format epsilon as '0.5' / '0.7' / '0.9' for filenames."""
    return f"{e:.2f}".rstrip("0").rstrip(".") if e != int(e) else f"{int(e)}"


def _fmt_d(d: float) -> str:
    """Format delta with 4 decimals, trimming trailing zeros (filename-safe)."""
    s = f"{d:.4f}"
    s = s.rstrip("0").rstrip(".")
    return s if s else "0"


def run_simplify(idv: int, delta: float, epsilon: float, baseline_name: str,
                  out_suffix: str = "", timeout=300, mem_cap_mb=0) -> Tuple[str, float]:
    """Run simplify and output to data/<id>/<out_name>.

    Returns (status, elapsed_seconds). status is one of 'ok', 'timeout',
    'memcap', 'error'. The elapsed time is measured around the headless
    subprocess with the monotonic clock.

    If `out_suffix` is empty, the file is named
    `simplify_against_<BASELINE>_<e>_<d>.txt` (per (e, d) invocation).
    Otherwise it is named `<out_suffix>` (e.g. legacy single-name form).
    """
    t0 = time.monotonic()
    data_dir = DATA_DIR / str(idv)
    if out_suffix:
        sh_out = data_dir / out_suffix
    else:
        sh_out = data_dir / f"simplify_against_{baseline_name}_{_fmt_eps(epsilon)}_{_fmt_d(delta)}.txt"
    cmd = [str(SIMPLIFY), str(idv), "--out", "-d", str(delta), "-e", str(epsilon)]
    status, _stdout, _stderr = run_cmd(cmd, cwd=REPO_ROOT, timeout=timeout,
                                        mem_cap_mb=mem_cap_mb, tag=f"simplify id={idv} d={delta:.2f} e={epsilon}")
    elapsed = time.monotonic() - t0
    if status != 'ok':
        # If the child was killed mid-run, the partial file it was writing
        # is still on disk and would be picked up by the next attempt as if
        # it were a valid output. Wipe it so a retry can produce a clean one.
        src = data_dir / "simplify.txt"
        if src.exists():
            try:
                src.unlink()
            except OSError:
                pass
        return status, elapsed
    src = data_dir / "simplify.txt"
    if src.exists():
        shutil.copy2(src, sh_out)
        return 'ok', elapsed
    return 'error', elapsed


def _median_run(cmd, *, cwd, timeout, mem_cap_mb, tag, n_repeats, on_success=None):
    """Run `cmd` up to n_repeats times in fresh processes, return median elapsed.

    Each repeat is a brand-new subprocess — there's no shared state between
    repeats, so the measurement amortizes dylib/Qt loaders over the first
    repeat and reports the typical run cost over the remainder. On Linux
    `fork()` would re-use a warm page cache, so we explicitly use
    `start_new_session=True`-style isolation by spawning fresh; on macOS
    POSIX `fork` is what `subprocess.Popen` uses by default and is cheap.

    on_success: optional callable invoked after each successful run with
    the elapsed seconds; use it to perform per-run post-processing
    (e.g. move a sidecar output file to its canonical location). If the
    callable raises, the run is still considered 'ok' for the purposes
    of timing — but the exception propagates.

    Returns (status, median_seconds) where median_seconds is the median
    of elapsed times across the repeats, or 0.0 if every repeat failed.
    """
    n = max(1, int(n_repeats))
    times = []
    last_status = 'ok'
    for i in range(n):
        t0 = time.monotonic()
        status, _, _ = run_cmd(cmd, cwd=cwd, timeout=timeout,
                                mem_cap_mb=mem_cap_mb, tag=f"{tag} rep{i+1}/{n}")
        elapsed = time.monotonic() - t0
        if status != 'ok':
            last_status = status
            # Don't pollute the median with a failed run.
            continue
        if on_success is not None:
            on_success()
        times.append(elapsed)
    if not times:
        return last_status, 0.0
    return 'ok', statistics.median(times)


def run_dp(original: Path, dp_out: Path, timeout=30, mem_cap_mb=0, n_repeats=1) -> Tuple[str, float]:
    cmd = [str(DP_TOOL), str(original), "200", str(dp_out)]
    return _median_run(cmd, cwd=REPO_ROOT, timeout=timeout, mem_cap_mb=mem_cap_mb,
                        tag=f"DP {original.parent.name}", n_repeats=n_repeats)


def run_dots(idv: int, dots_out: Path, timeout=60, mem_cap_mb=0, n_repeats=1) -> Tuple[str, float]:
    cmd = [str(DOTS_BIN), str(idv), "1000"]

    def _post():
        # DOTS writes to data/<id>/dots_simplified.txt (hardcoded by the C++ binary).
        # Move (not copy) it to the canonical DOTS.txt so we don't leave a
        # duplicate around — the dots_simplified.txt name leaks an internal
        # detail of the original DOTS C++ tool that callers shouldn't have to
        # know about.
        src = DATA_DIR / str(idv) / "dots_simplified.txt"
        if src.exists():
            shutil.move(str(src), str(dots_out))

    return _median_run(cmd, cwd=REPO_ROOT, timeout=timeout, mem_cap_mb=mem_cap_mb,
                        tag=f"DOTS id={idv}", n_repeats=n_repeats, on_success=_post)


def run_squish(original: Path, squish_out: Path, timeout=30, mem_cap_mb=0, n_repeats=1) -> Tuple[str, float]:
    cmd = [str(SQUISH_BIN), str(original), str(SQUISH_BUFFER), str(squish_out)]
    return _median_run(cmd, cwd=REPO_ROOT, timeout=timeout, mem_cap_mb=mem_cap_mb,
                        tag=f"SQUISH {original.parent.name}", n_repeats=n_repeats)


def needs_header(path: Path, expected_header: list[str]) -> bool:
    try:
        if not path.exists():
            return True
        with path.open('r', newline='') as f:
            first = f.readline()
            if not first:
                return True
            return first.strip() != ','.join(expected_header)
    except Exception:
        return True


def _id(row) -> Optional[int]:
    try:
        return int(row[0])
    except (ValueError, IndexError, TypeError):
        return None


def _algo(row) -> Optional[str]:
    try:
        return row[1]
    except (IndexError, TypeError):
        return None


def _load_existing_csv(path: Path):
    """Read the existing CSV and classify each (id, baseline_algo) row.

    Returns a tuple (valid_keys, reuse_data) where:
      * valid_keys:  set[(int, str)] of rows whose simplify_points and
                     simplify_frechet are both non-(-1). These rows are
                     considered complete and are NOT re-run on --resume.
      * reuse_data:  dict[(int, str), dict] of rows that ARE invalid (any -1
                     in the simplify_* fields). The dict carries the existing
                     baseline_points / baseline_frechet / baseline_time so
                     process_id can skip re-running the baseline algorithm
                     and jump straight to the simplify sweep.

    A row is also classified as invalid if it was emitted by an earlier
    baseline-skipped path (baseline_points == -1); the baseline is then
    re-run from scratch because we have no usable Frechet distance to derive
    d from.
    """
    valid_keys: set = set()
    reuse_data: dict = {}
    if not path.exists():
        return valid_keys, reuse_data
    try:
        with path.open('r', newline='') as f:
            reader = csv.reader(f)
            header = next(reader, None)
            if header is None:
                return valid_keys, reuse_data
            for row in reader:
                if not row or len(row) < 7:
                    continue
                try:
                    idv = int(row[0])
                    algo = row[1]
                    baseline_points = int(row[3])
                    simplify_points = int(row[4])
                    baseline_frechet = float(row[5])
                    simplify_frechet = float(row[6])
                except (ValueError, IndexError):
                    continue
                key = (idv, algo)
                baseline_ok = baseline_points != -1 and baseline_frechet != -1
                simplify_ok = simplify_points != -1 and simplify_frechet != -1
                if baseline_ok and simplify_ok:
                    valid_keys.add(key)
                elif baseline_ok:
                    # SH produced nothing or its Frechet failed, but the
                    # baseline side is intact — reuse baseline, rerun SH.
                    reuse_data[key] = {
                        'baseline_points': baseline_points,
                        'baseline_frechet': baseline_frechet,
                        'baseline_algo_time': float(row[7]) if len(row) > 7 and row[7] not in ('', '-1') else 0.0,
                    }
                # else: baseline itself failed (baseline_points == -1) —
                # nothing to reuse, leave the key out of both sets so
                # process_id runs the full pipeline for this (id, algo).
    except Exception as e:
        print(f"[resume] warning: failed to read {path}: {e}", file=sys.stderr)
    return valid_keys, reuse_data


def _csv_row_for_key(path: Path, key: tuple) -> Optional[dict]:
    """Return a reuse_data-style dict for the row matching (id, algo) in `path`.

    Returns None if the row is missing or unparseable. Used by
    --recompute-baseline to convert a "valid" row into a "reuse" row
    (so process_id re-times the baseline while keeping the simplify
    sweep's output).
    """
    if not path.exists():
        return None
    try:
        with path.open('r', newline='') as f:
            reader = csv.reader(f)
            next(reader, None)  # header
            for row in reader:
                if not row or len(row) < 7:
                    continue
                try:
                    if (int(row[0]), row[1]) != key:
                        continue
                    baseline_points = int(row[3])
                    baseline_frechet = float(row[5])
                except (ValueError, IndexError):
                    continue
                if baseline_points == -1 or baseline_frechet == -1:
                    return None
                return {
                    'baseline_points': baseline_points,
                    'baseline_frechet': baseline_frechet,
                    'baseline_algo_time': float(row[7]) if len(row) > 7 and row[7] not in ('', '-1') else 0.0,
                }
    except Exception as e:
        print(f"[resume] warning: failed to read {path}: {e}", file=sys.stderr)
    return None


def _remove_stale_rows(path: Path, keys_to_remove: set) -> None:
    """Rewrite `path` dropping rows whose (id, algo) is in keys_to_remove.

    Used when a previously-invalid (id, algo) row is being replaced with a
    fresh run — we drop the old -1 row so the CSV doesn't end up with
    duplicates. Caller must hold the write lock (the main process is the
    only writer, so this is naturally serialized).
    """
    if not keys_to_remove or not path.exists():
        return
    try:
        with path.open('r', newline='') as f:
            reader = csv.reader(f)
            header = next(reader, None)
            kept = [r for r in reader if r and (int(r[0]), r[1]) not in keys_to_remove]
    except Exception as e:
        print(f"[resume] warning: failed to scan {path} for cleanup: {e}",
              file=sys.stderr)
        return
    tmp = path.with_suffix(path.suffix + '.tmp')
    try:
        with tmp.open('w', newline='') as f:
            writer = csv.writer(f)
            if header is not None:
                writer.writerow(header)
            writer.writerows(kept)
        tmp.replace(path)
    except Exception as e:
        print(f"[resume] warning: failed to rewrite {path}: {e}",
              file=sys.stderr)
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass


def process_id(idv: int, mem_cap_mb: int = 0,
                algo_timeout: int = DEFAULT_ALGO_TIMEOUT_S,
                sh_timeout: int = DEFAULT_SH_TIMEOUT_S,
                frechet_timeout: int = DEFAULT_FRECHET_TIMEOUT_S,
                n_repeats: int = 1,
                valid_keys: Optional[set] = None,
                reuse_data: Optional[dict] = None) -> list:
    """Process one ID and return rows for CSV.

    valid_keys / reuse_data are populated by main() from the existing
    compare_points.csv on --resume. For each (id, baseline_algo):
      * if the key is in valid_keys: skip — row is already complete.
      * elif the key is in reuse_data: skip the baseline, reuse the
        existing baseline_points / baseline_frechet / baseline_algo_time,
        and re-run the simplify sweep only.
      * else: full pipeline (baseline + simplify).
    Rows whose (id, algo) is reused are emitted with a leading 'R' marker
    in the row's last column so downstream tools can tell them apart, but
    the schema is unchanged.

    n_repeats controls how many times each baseline binary is run per
    (id, baseline); the recorded baseline_time is the median of the
    per-invocation wall times. (default 1; --repeats N on the CLI sets
    it to 3 by default.) simplify_time is the subprocess wall time for
    one headless simplify invocation, so it does not depend on n_repeats.
    """
    data_dir = DATA_DIR / str(idv)
    original = data_dir / "original.txt"

    if not original.exists():
        print(f"[ID:{idv}] Missing original.txt, skipping", file=sys.stderr)
        return []

    orig_points = count_points(original)
    print(f"[ID:{idv}] Starting...")
    rows = []

    baselines = {
        'DP': {
            'file': data_dir / "DP.txt",
            'run': lambda: run_dp(original, data_dir / "DP.txt",
                                  timeout=algo_timeout, mem_cap_mb=mem_cap_mb,
                                  n_repeats=n_repeats),
        },
        'DOTS': {
            'file': data_dir / "DOTS.txt",
            'run': lambda: run_dots(idv, data_dir / "DOTS.txt",
                                    timeout=algo_timeout, mem_cap_mb=mem_cap_mb,
                                    n_repeats=n_repeats),
        },
        'SQUISH': {
            'file': data_dir / "SQUISH.txt",
            'run': lambda: run_squish(original, data_dir / "SQUISH.txt",
                                      timeout=algo_timeout, mem_cap_mb=mem_cap_mb,
                                      n_repeats=n_repeats),
        },
    }

    valid_keys = valid_keys or set()
    reuse_data = reuse_data or {}

    for algo_name, algo in baselines.items():
        key = (idv, algo_name)

        # Resume: skip a row that is already complete.
        if key in valid_keys:
            print(f"[ID:{idv}] {algo_name} already valid in CSV, skipping")
            continue

        # Resume: reuse the existing baseline result and re-run the
        # simplify sweep only. The baseline algorithm is not invoked
        # again — we trust the previously-recorded baseline_points /
        # baseline_frechet / baseline_algo_time and derive a fresh delta
        # from the existing Frechet distance.
        if key in reuse_data:
            cached = reuse_data[key]
            baseline_points = cached['baseline_points']
            baseline_fret = cached['baseline_frechet']
            baseline_algo_time = cached['baseline_algo_time']
            baseline_fret_time = 0.0
            print(f"[ID:{idv}] {algo_name} reusing cached baseline "
                  f"(points={baseline_points}, frechet={baseline_fret}, "
                  f"time={baseline_algo_time:.2f}s); re-running simplify")
            rows.extend(_sweep_simplify(
                idv, algo_name, orig_points, baseline_points, baseline_fret,
                baseline_algo_time, baseline_fret_time,
                sh_timeout=sh_timeout, frechet_timeout=frechet_timeout,
                mem_cap_mb=mem_cap_mb, baseline_algo_file=None))
            continue

        print(f"[ID:{idv}] Running {algo_name}...")
        status, baseline_algo_time = algo['run']()
        algo_file = algo['file']

        if status != 'ok' or not algo_file.exists():
            print(f"[ID:{idv}]   {algo_name} {status}, skipping")
            rows.append([idv, algo_name, orig_points, -1, -1, -1, -1, baseline_algo_time, -1, -1, -1])
            continue

        baseline_points = count_points(algo_file)
        print(f"[ID:{idv}]   {algo_name} done, points={baseline_points}, time={baseline_algo_time:.2f}s")
        print(f"[ID:{idv}]   Computing Frechet for {algo_name}...")
        baseline_fret, baseline_fret_time, fret_status = frechet(
            original, algo_file, timeout=frechet_timeout, mem_cap_mb=mem_cap_mb,
            tag=f"frechet id={idv} {algo_name}")
        print(f"[ID:{idv}]   {algo_name} Frechet = {baseline_fret} (time={baseline_fret_time:.2f}s, status={fret_status})")

        if baseline_fret is None:
            print(f"[ID:{idv}]   {algo_name} Frechet {fret_status}, skipping")
            rows.append([idv, algo_name, orig_points, baseline_points, -1, baseline_fret, -1,
                        baseline_algo_time, -1, -1, -1])
            continue

        # Batch Frechet for simplify: run all 3 simplify first, then one Julia call
        # (baseline Frechet must come first since it's used to derive delta)
        simp_results = []
        for eps in EPS_VALUES:
            delta = baseline_fret / (1.0 + eps)
            print(f"[ID:{idv}]   Running simplify e={eps} d={delta:.4f}...")
            sh_status, simplify_time = run_simplify(
                idv, delta, eps, algo_name,
                timeout=sh_timeout, mem_cap_mb=mem_cap_mb)
            sh_file = DATA_DIR / str(idv) / f"simplify_against_{algo_name}_{_fmt_eps(eps)}_{_fmt_d(delta)}.txt"
            if sh_status != 'ok' or not sh_file.exists():
                print(f"[ID:{idv}]     simplify e={eps} {sh_status}, skipping")
                simp_results.append((eps, delta, None, sh_status, simplify_time, None))
            else:
                simplify_points = count_points(sh_file)
                print(f"[ID:{idv}]     simplify e={eps} done, points={simplify_points}, time={simplify_time:.2f}s")
                simp_results.append((eps, delta, sh_file, 'ok', simplify_time, simplify_points))

        # Batch all simplify Frechet calls into ONE Julia invocation
        batch_paths = [sh_file for _, _, sh_file, status, _, _ in simp_results
                       if status == 'ok' and sh_file is not None]

        if not batch_paths:
            rows.append([idv, algo_name, orig_points, baseline_points, -1, baseline_fret, -1,
                        baseline_algo_time, -1, -1, -1])
            continue

        print(f"[ID:{idv}]   Batch Frechet for {len(batch_paths)} simplify paths...")
        batch_results, batch_total_time, batch_status = frechet_batch(
            idv, batch_paths, timeout=frechet_timeout, mem_cap_mb=mem_cap_mb,
            tag=f"frechet_batch id={idv} {algo_name}")
        fret_map = {str(p): (d, s) for p, d, s in batch_results}

        # Pick best simplify result
        best_simp = None
        best_simp_fret = None
        best_simp_time = None
        best_e = None
        best_d = None

        for eps, delta, sh_file, status, simplify_time, simplify_points in simp_results:
            if status != 'ok' or sh_file is None:
                print(f"[ID:{idv}]     simplify e={eps} Frechet = (skipped)")
                continue
            sh_fret, sh_fret_status = fret_map.get(str(sh_file), (None, 'missing'))
            print(f"[ID:{idv}]     simplify e={eps} Frechet = {sh_fret} (batch, status={sh_fret_status})")
            if sh_fret is not None and (best_simp_fret is None or sh_fret < best_simp_fret):
                best_simp_fret = sh_fret
                best_simp = simplify_points
                # simplify_time is the headless subprocess wall time. The
                # Frechet cost is excluded from this column; that call only
                # selects the best (e, d).
                best_simp_time = simplify_time
                best_e = eps
                best_d = delta

        if best_simp is None:
            rows.append([idv, algo_name, orig_points, baseline_points, -1, baseline_fret, -1,
                        baseline_algo_time, -1, -1, -1])
        else:
            rows.append([
                idv, algo_name,
                orig_points,
                baseline_points, best_simp,
                baseline_fret, best_simp_fret,
                baseline_algo_time, best_simp_time,
                best_e, best_d,
            ])

    print(f"[ID:{idv}] Done.")
    return rows


def _sweep_simplify(idv: int, algo_name: str, orig_points, baseline_points,
                       baseline_fret, baseline_algo_time, baseline_fret_time,
                       sh_timeout, frechet_timeout, mem_cap_mb,
                       baseline_algo_file: Optional[Path] = None) -> list:
    """Run the simplify 3-epsilon sweep for one (id, baseline) and return
    one CSV row (or one all-(-1) row if every simplify attempt failed).

    Runs all 3 simplify calls first, then batches all 4 Frechet calls
    (baseline + 3 simplify) into a single Julia invocation.
    """
    # Phase 1: run all 3 simplify invocations, collecting results
    simp_results = []  # [(eps, delta, sh_file, sh_status, simplify_time, simplify_points), ...]

    for eps in EPS_VALUES:
        delta = baseline_fret / (1.0 + eps)
        print(f"[ID:{idv}]   Running simplify e={eps} d={delta:.4f}...")
        sh_status, simplify_time = run_simplify(
            idv, delta, eps, algo_name,
            timeout=sh_timeout, mem_cap_mb=mem_cap_mb)
        sh_file = DATA_DIR / str(idv) / f"simplify_against_{algo_name}_{_fmt_eps(eps)}_{_fmt_d(delta)}.txt"

        if sh_status != 'ok' or not sh_file.exists():
            print(f"[ID:{idv}]     simplify e={eps} {sh_status}, skipping")
            simp_results.append((eps, delta, None, sh_status, simplify_time, None))
        else:
            simplify_points = count_points(sh_file)
            print(f"[ID:{idv}]     simplify e={eps} done, points={simplify_points}, time={simplify_time:.2f}s")
            simp_results.append((eps, delta, sh_file, 'ok', simplify_time, simplify_points))

    # Phase 2: batch all Frechet calls into ONE Julia invocation
    # Build paths list: baseline_algo_file + all successful simplify files
    batch_paths = []
    if baseline_algo_file is not None and baseline_algo_file.exists():
        batch_paths.append(baseline_algo_file)
    for _, _, sh_file, status, _, _ in simp_results:
        if status == 'ok' and sh_file is not None:
            batch_paths.append(sh_file)

    if not batch_paths:
        return [[idv, algo_name, orig_points, baseline_points, -1, baseline_fret, -1,
                baseline_algo_time, -1, -1, -1]]

    print(f"[ID:{idv}]   Batch Frechet for {len(batch_paths)} paths...")
    batch_results, batch_total_time, batch_status = frechet_batch(
        idv, batch_paths, timeout=frechet_timeout, mem_cap_mb=mem_cap_mb,
        tag=f"frechet_batch id={idv} {algo_name}")

    # Parse batch results: map path -> (distance, status)
    fret_map = {}
    for path, dist, st in batch_results:
        fret_map[path] = (dist, st)

    # Extract baseline Frechet (first in batch)
    baseline_fret_new = None
    baseline_fret_time_new = 0.0
    if baseline_algo_file and str(baseline_algo_file) in fret_map:
        baseline_fret_new, bst = fret_map[str(baseline_algo_file)]
        baseline_fret_time_new = batch_total_time / len(batch_paths) if batch_status == 'ok' else 0.0
        print(f"[ID:{idv}]   {algo_name} Frechet = {baseline_fret_new} (batch, status={bst})")

    # Extract simplify Frechet results and find best
    best_simp = None
    best_simp_fret = None
    best_simp_time = None
    best_e = None
    best_d = None

    for eps, delta, sh_file, status, simplify_time, simplify_points in simp_results:
        if status != 'ok' or sh_file is None:
            print(f"[ID:{idv}]     simplify e={eps} Frechet = (baseline only)")
            continue

        sh_fret, sh_fret_status = fret_map.get(str(sh_file), (None, 'missing'))
        print(f"[ID:{idv}]     simplify e={eps} Frechet = {sh_fret} (batch, status={sh_fret_status})")

        if sh_fret is not None and (best_simp_fret is None or sh_fret < best_simp_fret):
            best_simp_fret = sh_fret
            best_simp = simplify_points
            # Headless simplify subprocess wall time. Frechet is intentionally
            # excluded from the simplify_time column (see process_id).
            best_simp_time = simplify_time
            best_e = eps
            best_d = delta

    if best_simp is None:
        return [[idv, algo_name, orig_points, baseline_points, -1, baseline_fret, -1,
                baseline_algo_time, -1, -1, -1]]
    return [[
        idv, algo_name,
        orig_points,
        baseline_points, best_simp,
        baseline_fret if baseline_fret_new is None else baseline_fret_new,
        best_simp_fret,
        baseline_algo_time, best_simp_time,
        best_e, best_d,
    ]]


def _process_id_worker(idv, mem_cap_mb, algo_timeout, sh_timeout,
                       frechet_timeout, n_repeats, valid_keys, reuse_data):
    """Top-level worker for ProcessPoolExecutor.

    Must be module-level (not nested in main) so it can be pickled by the
    default multiprocessing pickler. functools.partial binds the
    per-run arguments and is itself picklable.
    """
    return process_id(
        idv,
        mem_cap_mb=mem_cap_mb,
        algo_timeout=algo_timeout,
        sh_timeout=sh_timeout,
        frechet_timeout=frechet_timeout,
        n_repeats=n_repeats,
        valid_keys=valid_keys,
        reuse_data=reuse_data,
    )


def main():
    parser = argparse.ArgumentParser(description="Compare simplify against DP, DOTS, SQUISH")
    parser.add_argument('--a', type=int, default=1, help='start id (inclusive)')
    parser.add_argument('--b', type=int, default=1000, help='end id (inclusive)')
    parser.add_argument('--resume', action='store_true',
                        help='read existing compare_points.csv and skip (id, baseline_algo) rows '
                             'that are already valid (no -1 in simplify_points/simplify_frechet). '
                             'Rows with -1 in those columns are re-run, reusing the cached baseline '
                             'numbers and only invoking simplify again. Without this flag the '
                             'benchmark always re-runs from scratch.')
    parser.add_argument('--workers', '-w', type=int, default=8,
                        help='number of parallel worker processes (default: 8)')
    parser.add_argument('--mem-cap-mb', type=int, default=DEFAULT_MEM_CAP_MB,
                        help=f'RSS cap per child process in MB (default: {DEFAULT_MEM_CAP_MB}; '
                             '0 disables the watchdog)')
    parser.add_argument('--algo-timeout', type=int, default=DEFAULT_ALGO_TIMEOUT_S,
                        help=f'wall-clock timeout in seconds for DP/DOTS/SQUISH '
                             f'(default: {DEFAULT_ALGO_TIMEOUT_S})')
    parser.add_argument('--sh-timeout', type=int, default=DEFAULT_SH_TIMEOUT_S,
                        help=f'wall-clock timeout in seconds per simplify invocation '
                             f'(default: {DEFAULT_SH_TIMEOUT_S})')
    parser.add_argument('--frechet-timeout', type=int, default=DEFAULT_FRECHET_TIMEOUT_S,
                        help=f'wall-clock timeout in seconds per Frechet (Julia) call '
                             f'(default: {DEFAULT_FRECHET_TIMEOUT_S})')
    parser.add_argument('--repeats', type=int, default=3,
                        help='number of times to run each baseline binary (DP / DOTS / SQUISH) '
                             'per (id, baseline). The recorded baseline_time is the median of '
                             'the per-invocation wall times; the first repeat amortizes loader '
                             'cost. simplify_time is one headless subprocess wall time and is '
                             'independent of this flag. (default: 3)')
    parser.add_argument('--recompute-baseline', action='store_true',
                        help='on --resume, re-run the baseline binary for every (id, baseline) '
                             'even if the row is fully valid. Keeps the cached simplify result '
                             'and Frechet calls but refreshes the baseline_time column. Use this '
                             'after changing --repeats to migrate the existing CSV to the new '
                             'timing policy.')
    args = parser.parse_args()

    if args.repeats < 1:
        parser.error("--repeats must be >= 1")

    # Get IDs to process
    ids = list(range(args.a, args.b + 1))

    # On --resume, load existing CSV state and decide per (id, algo) what
    # needs work. The old behavior treated the whole ID as completed if any
    # row for it existed in the CSV — which meant a row with -1 in the
    # simplify columns would never be retried. The new behavior is row-
    # level: a row is "valid" only when both simplify_points and
    # simplify_frechet are real numbers (not -1), and an ID is skipped
    # entirely only if all 3 of its rows are valid.
    valid_keys: set = set()
    reuse_data: dict = {}
    if args.resume and OUTPUT_CSV.exists():
        valid_keys, reuse_data = _load_existing_csv(OUTPUT_CSV)
        if args.recompute_baseline:
            # Move every currently-valid key into "reuse" so the baseline
            # is re-timed while the simplify sweep result is preserved.
            # We need the cached row to be available to _sweep_simplify,
            # so copy it from valid_keys first by reading the CSV.
            for key in list(valid_keys):
                reuse_data[key] = _csv_row_for_key(OUTPUT_CSV, key)
            valid_keys = set()
            print(f"Recompute baseline: forced {len(reuse_data)} (id, algo) rows to re-time the baseline")
        skipped_rows = len(valid_keys)
        invalid_rows = len(reuse_data)
        # An ID is "needs work" if at least one of its 3 algorithm rows is
        # missing or invalid. We map over ids in input order (not after
        # filtering) so the count below matches what users see in the log.
        needs_work = lambda idv: any(
            (idv, a) not in valid_keys for a in ('DP', 'DOTS', 'SQUISH')
        )
        skipped_ids = sum(1 for i in ids if not needs_work(i))
        ids = [i for i in ids if needs_work(i)]
        print(f"Resuming: {skipped_rows} valid rows across {skipped_ids} fully-done IDs skipped; "
              f"{invalid_rows} invalid rows queued for retry; {len(ids)} IDs to process")

    cols = ['id', 'baseline_algo', 'original_points', 'baseline_points', 'simplify_points',
            'baseline_frechet', 'simplify_frechet',
            'baseline_time', 'simplify_time', 'best_e', 'best_d']

    # Write header only if file doesn't exist
    write_header = not OUTPUT_CSV.exists()
    if write_header:
        with OUTPUT_CSV.open('w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(cols)

    print(f"Processing IDs {args.a}-{args.b} ({len(ids)} IDs) with {args.workers} workers...")
    print(f"Resource limits: mem_cap={args.mem_cap_mb} MB, "
          f"algo_timeout={args.algo_timeout}s, sh_timeout={args.sh_timeout}s, "
          f"frechet_timeout={args.frechet_timeout}s, repeats={args.repeats}")

    _process = functools.partial(
        _process_id_worker,
        mem_cap_mb=args.mem_cap_mb,
        algo_timeout=args.algo_timeout,
        sh_timeout=args.sh_timeout,
        frechet_timeout=args.frechet_timeout,
        n_repeats=args.repeats,
        valid_keys=valid_keys,
        reuse_data=reuse_data,
    )

    def _writeback_rows(rows: list) -> None:
        """Append rows to the CSV, dropping any prior row for the same
        (id, algo) so the CSV stays at one row per key (the latest one
        wins). Caller must be the only writer (the main process).
        """
        if not rows:
            return
        touched = {(_id(r), _algo(r)) for r in rows if _id(r) is not None}
        if touched:
            _remove_stale_rows(OUTPUT_CSV, touched)
        with OUTPUT_CSV.open('a', newline='') as f:
            writer = csv.writer(f)
            for row in rows:
                writer.writerow(row)

    total_rows = 0
    completed = 0
    if args.workers == 1:
        # Serial mode — useful when ProcessPoolExecutor hits sandbox/permission
        # issues (e.g. macOS sysconf SC_SEM_NSEMS_MAX denied).
        for idv in ids:
            try:
                rows = _process(idv)
            except Exception as e:
                print(f"[ID:{idv}] worker raised: {e}", file=sys.stderr)
                rows = []
            _writeback_rows(rows)
            total_rows += len(rows)
            completed += 1
            if completed % 1 == 0 or completed == len(ids):
                print(f"Progress: {completed}/{len(ids)} IDs processed ({total_rows} rows)")
        return 0

    # Use fork() on macOS to avoid the spawn-based pickling failure that
    # happens when functools.partial of a top-level worker is submitted to
    # a default ProcessPoolExecutor. The child doesn't hold any C++/Qt
    # state — workers only subprocess out to external binaries — so fork
    # is safe here.
    import multiprocessing as _mp
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=args.workers, mp_context=_mp.get_context('fork')
    ) as executor:
        # Use as_completed so we get progress feedback as IDs finish,
        # regardless of completion order. Each worker writes to its own
        # data/<id>/ dir, so the per-ID writes are safe to do inside
        # process_id(); we just collect the resulting rows and append to
        # the CSV here in the main process.
        future_to_id = {executor.submit(_process, idv): idv for idv in ids}
        for future in concurrent.futures.as_completed(future_to_id):
            idv = future_to_id[future]
            try:
                rows = future.result()
            except Exception as e:
                print(f"[ID:{idv}] worker raised: {e}", file=sys.stderr)
                rows = []
            _writeback_rows(rows)
            total_rows += len(rows)
            completed += 1
            if completed % 10 == 0 or completed == len(ids):
                print(f"Progress: {completed}/{len(ids)} IDs processed ({total_rows} rows)")

    print(f"Done. {total_rows} rows written to {OUTPUT_CSV}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
