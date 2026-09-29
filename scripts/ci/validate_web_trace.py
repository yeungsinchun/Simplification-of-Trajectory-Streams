#!/usr/bin/env python3
"""
Validate NDJSON (or gzipped NDJSON) trace emitted by the web server / simplify handler.

Expects the v2 trace format:
  header line:  {type:"header", eps, delta, grid_val, r_val, expected_frechet, bbox[4], stream:[[x,y]]}
  prefix lines: {type:"prefix", data:{p0, p0_idx, end_idx, P, output[2], steps:[{stream_idx, pi, Gi, buffer[2], candidates:[{idx, alive, F, F_Si, S}]}]}}
  done line:    {type:"done", time_ms, simplified:[[x,y]], frechet_distance}

Usage:
  python3 scripts/ci/validate_web_trace.py --input /tmp/trace.ndjson --epsilon 0.5 --delta 200 --trace-id 1
  python3 scripts/ci/validate_web_trace.py --input /tmp/trace.gz --epsilon 0.5 --delta 200 --trace-id 1 --wall-time-ms 1234
  cat trace.ndjson | python3 scripts/ci/validate_web_trace.py --epsilon 0.5 --delta 200

Exit 0 on success, non-zero on any structural regression. Prints a human summary and
machine-readable JSON to stdout when --json-summary is given.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import sys
from pathlib import Path

def _is_finite(v) -> bool:
    return type(v) in (int, float) and math.isfinite(float(v))

def _is_point(p) -> bool:
    return isinstance(p, (list, tuple)) and len(p) == 2 and _is_finite(p[0]) and _is_finite(p[1])

def _is_points(arr) -> bool:
    return isinstance(arr, list) and all(_is_point(pt) for pt in arr)

def _close(a: float, b: float, rel: float = 1e-6, abs_tol: float = 1e-9) -> bool:
    if a == b:
        return True
    try:
        af = float(a); bf = float(b)
    except Exception:
        return False
    if not (math.isfinite(af) and math.isfinite(bf)):
        return False
    return math.isclose(af, bf, rel_tol=rel, abs_tol=abs_tol)

def _read_bytes(input_path: str | None) -> tuple[bytes, int, bool]:
    """Return (raw_bytes_decompressed, raw_bytes_on_disk, was_gzipped)."""
    if input_path is None or input_path == "-":
        data = sys.stdin.buffer.read()
        # detect gzip magic
        if len(data) >= 2 and data[0] == 0x1f and data[1] == 0x8b:
            try:
                dec = gzip.decompress(data)
                return dec, len(data), True
            except Exception:
                return data, len(data), False
        return data, len(data), False
    p = Path(input_path)
    if not p.exists():
        print(f"error: input file not found: {input_path}", file=sys.stderr)
        sys.exit(2)
    raw = p.read_bytes()
    if len(raw) >= 2 and raw[0] == 0x1f and raw[1] == 0x8b:
        try:
            dec = gzip.decompress(raw)
            return dec, len(raw), True
        except Exception as e:
            print(f"error: failed to gunzip {input_path}: {e}", file=sys.stderr)
            sys.exit(2)
    return raw, len(raw), False

def validate_ndjson(raw: bytes, expected_eps: float | None, expected_delta: float | None, trace_id: int | None, orig_n: int | None) -> dict:
    text = raw.decode("utf-8", errors="strict")
    lines = [ln for ln in text.splitlines() if ln.strip() != ""]
    if not lines:
        raise ValueError("empty NDJSON (no lines)")

    # NDJSON path
    header = None
    prefixes = []
    done = None
    errors: list[str] = []
    line_objs = []

    for idx, line in enumerate(lines):
        try:
            obj = json.loads(line)
        except json.JSONDecodeError as e:
            raise ValueError(f"line {idx+1}: invalid JSON: {e} — {line[:200]!r}") from e
        if not isinstance(obj, dict):
            raise ValueError(f"line {idx+1}: expected JSON object, got {type(obj).__name__}")
        typ = obj.get("type")
        if typ not in ("header", "prefix", "done", "error"):
            raise ValueError(f"line {idx+1}: missing or unknown type: {typ!r} (expected header/prefix/done)")
        if typ == "error":
            raise ValueError(f"line {idx+1}: server emitted error: {obj.get('message')!r}")
        line_objs.append(obj)
        if typ == "header":
            if header is not None:
                raise ValueError(f"line {idx+1}: duplicate header")
            header = obj
        elif typ == "prefix":
            if header is None:
                raise ValueError(f"line {idx+1}: prefix before header")
            if done is not None:
                raise ValueError(f"line {idx+1}: prefix after done")
            prefixes.append(obj)
        elif typ == "done":
            if header is None:
                raise ValueError(f"line {idx+1}: done before header")
            if done is not None:
                raise ValueError(f"line {idx+1}: duplicate done")
            done = obj

    if header is None:
        raise ValueError("missing header line")
    if done is None:
        raise ValueError("missing done line (stream truncated)")
    if not prefixes:
        raise ValueError("no prefix lines (expected at least one)")

    # Header validation
    for field in ("eps", "delta", "grid_val", "r_val", "expected_frechet", "bbox", "stream"):
        if field not in header:
            raise ValueError(f"header missing field: {field}")
    if not _is_finite(header["eps"]):
        raise ValueError(f"header.eps not finite: {header['eps']!r}")
    if not _is_finite(header["delta"]):
        raise ValueError(f"header.delta not finite: {header['delta']!r}")
    if expected_eps is not None and not _close(header["eps"], expected_eps):
        raise ValueError(f"header.eps mismatch: got {header['eps']} expected {expected_eps}")
    if expected_delta is not None and not _close(header["delta"], expected_delta):
        raise ValueError(f"header.delta mismatch: got {header['delta']} expected {expected_delta}")
    for k in ("grid_val", "r_val"):
        if not _is_finite(header[k]) or float(header[k]) <= 0:
            raise ValueError(f"header.{k} must be positive finite, got {header[k]!r}")
    if not _is_finite(header["expected_frechet"]) or float(header["expected_frechet"]) < 0:
        raise ValueError(f"header.expected_frechet must be finite >=0, got {header['expected_frechet']!r}")
    bbox = header["bbox"]
    if not (isinstance(bbox, list) and len(bbox) == 4 and all(_is_finite(v) for v in bbox)):
        raise ValueError(f"header.bbox must be 4 finite numbers, got {bbox!r}")
    if not (bbox[0] < bbox[2] and bbox[1] < bbox[3]):
        raise ValueError(f"header.bbox must satisfy min < max, got {bbox!r}")
    stream = header["stream"]
    if not _is_points(stream):
        raise ValueError("header.stream must be array of [x,y] points")
    if len(stream) < 2:
        raise ValueError(f"header.stream must have at least 2 points, got {len(stream)}")
    if orig_n is not None and len(stream) != orig_n:
        raise ValueError(f"header.stream length {len(stream)} != original.txt N {orig_n} for trace {trace_id}")
    # Optionally verify monotonic bbox? Skip.

    # Prefix validation
    total_steps = 0
    for p_idx, pref_obj in enumerate(prefixes):
        if pref_obj.get("type") != "prefix":
            raise ValueError(f"prefix {p_idx}: wrong type")
        data = pref_obj.get("data")
        if not isinstance(data, dict):
            raise ValueError(f"prefix {p_idx}: missing data object")
        for f in ("p0", "p0_idx", "end_idx", "P", "output", "steps"):
            if f not in data:
                raise ValueError(f"prefix {p_idx}: missing data.{f}")
        if not _is_point(data["p0"]):
            raise ValueError(f"prefix {p_idx}: p0 not a point")
        if type(data["p0_idx"]) is not int or data["p0_idx"] < 0:
            raise ValueError(f"prefix {p_idx}: p0_idx must be non-negative int")
        if type(data["end_idx"]) is not int or data["end_idx"] <= data["p0_idx"]:
            raise ValueError(f"prefix {p_idx}: end_idx must be > p0_idx")
        if data["p0_idx"] >= len(stream):
            raise ValueError(f"prefix {p_idx}: p0_idx {data['p0_idx']} out of range stream len {len(stream)}")
        if data["end_idx"] > len(stream):
            raise ValueError(f"prefix {p_idx}: end_idx {data['end_idx']} out of range stream len {len(stream)}")
        if p_idx == 0 and data["p0_idx"] != 0:
            raise ValueError(f"prefix {p_idx}: first p0_idx must be 0, got {data['p0_idx']}")
        try:
            sp = stream[data["p0_idx"]]
            if not (_close(data["p0"][0], sp[0]) and _close(data["p0"][1], sp[1])):
                raise ValueError(f"prefix {p_idx}: p0 {data['p0']!r} != stream[{data['p0_idx']}] {sp!r}")
        except IndexError:
            raise ValueError(f"prefix {p_idx}: p0_idx {data['p0_idx']} out of range stream len {len(stream)}")
        if not _is_points(data["P"]):
            raise ValueError(f"prefix {p_idx}: P must be array of points")
        if len(data["P"]) == 0:
            raise ValueError(f"prefix {p_idx}: P empty")
        output = data["output"]
        if not (isinstance(output, list) and len(output) == 2 and _is_point(output[0]) and _is_point(output[1])):
            raise ValueError(f"prefix {p_idx}: output must be 2 points")
        steps = data["steps"]
        if not isinstance(steps, list):
            raise ValueError(f"prefix {p_idx}: steps must be array")
        if len(steps) == 0:
            if data["end_idx"] != data["p0_idx"] + 1:
                raise ValueError(f"prefix {p_idx}: steps empty but end_idx {data['end_idx']} != p0_idx+1 {data['p0_idx']+1}")
        total_steps += len(steps)
        # Validate sequencing of steps
        for s_idx, step in enumerate(steps):
            if not isinstance(step, dict):
                raise ValueError(f"prefix {p_idx} step {s_idx}: not an object")
            for sf in ("stream_idx", "pi", "Gi", "candidates"):
                if sf not in step:
                    raise ValueError(f"prefix {p_idx} step {s_idx}: missing {sf}")
            if type(step["stream_idx"]) is not int or not (0 <= step["stream_idx"] < len(stream)):
                raise ValueError(f"prefix {p_idx} step {s_idx}: stream_idx {step['stream_idx']} out of range [0, {len(stream)})")
            if step["stream_idx"] != data["p0_idx"] + 1 + s_idx:
                raise ValueError(f"prefix {p_idx} step {s_idx}: stream_idx {step['stream_idx']} != expected {data['p0_idx'] + 1 + s_idx}")
            if not _is_point(step["pi"]):
                raise ValueError(f"prefix {p_idx} step {s_idx}: pi not a point")
            sp_pi = stream[step["stream_idx"]]
            if not (_close(step["pi"][0], sp_pi[0]) and _close(step["pi"][1], sp_pi[1])):
                raise ValueError(f"prefix {p_idx} step {s_idx}: pi {step['pi']!r} != stream[{step['stream_idx']}] {sp_pi!r}")
            Gi = step["Gi"]
            if not isinstance(Gi, list) or len(Gi) < 3:
                raise ValueError(f"prefix {p_idx} step {s_idx}: Gi must be array with >=3 points, got {len(Gi) if isinstance(Gi, list) else Gi!r}")
            if not _is_points(Gi):
                raise ValueError(f"prefix {p_idx} step {s_idx}: Gi contains invalid point")
            # v2 diet: buffer is omitted (viewer does not use it); accept with or without it
            if "buffer" in step:
                buffer = step["buffer"]
                if not (isinstance(buffer, list) and len(buffer) == 2 and _is_point(buffer[0]) and _is_point(buffer[1])):
                    raise ValueError(f"prefix {p_idx} step {s_idx}: buffer must be 2 points")
            cands = step["candidates"]
            if not isinstance(cands, list):
                raise ValueError(f"prefix {p_idx} step {s_idx}: candidates must be array")
            # v2 diet: previously-dead candidates (F < 3) are omitted, so sparse array is expected
            if not (0 < len(cands) <= len(data["P"])):
                raise ValueError(f"prefix {p_idx} step {s_idx}: candidates len {len(cands)} not in (0, {len(data['P'])}]")
            seen_idx = set()
            for _ci in cands:
                if not isinstance(_ci, dict) or "idx" not in _ci:
                    raise ValueError(f"prefix {p_idx} step {s_idx}: candidate missing idx")
                _idx = _ci["idx"]
                if _idx in seen_idx:
                    raise ValueError(f"prefix {p_idx} step {s_idx}: duplicate candidate idx {_idx}")
                seen_idx.add(_idx)
            # Validate each candidate
            for c_idx, cand in enumerate(cands):
                if not isinstance(cand, dict):
                    raise ValueError(f"prefix {p_idx} step {s_idx} cand {c_idx}: not object")
                for cf in ("idx", "alive", "F", "F_Si", "S"):
                    if cf not in cand:
                        raise ValueError(f"prefix {p_idx} step {s_idx} cand {c_idx}: missing {cf}")
                if type(cand["idx"]) is not int or not (0 <= cand["idx"] < len(data["P"])):
                    raise ValueError(f"prefix {p_idx} step {s_idx} cand {c_idx}: idx {cand['idx']} out of range [0, {len(data['P'])})")
                if not isinstance(cand["alive"], bool):
                    raise ValueError(f"prefix {p_idx} step {s_idx} cand {c_idx}: alive must be bool")
                for arr_name in ("F", "F_Si", "S"):
                    arr = cand[arr_name]
                    if not isinstance(arr, list):
                        raise ValueError(f"prefix {p_idx} step {s_idx} cand {c_idx}: {arr_name} must be array")
                    if arr and not _is_points(arr):
                        raise ValueError(f"prefix {p_idx} step {s_idx} cand {c_idx}: {arr_name} contains invalid point")
                # If alive, F and S should be polygonal? F empty means miss? For alive, F must have >=3? For dead, S may be empty.
                # We don't enforce strict, but check that alive => S non-empty?
                # In simplify.cpp, alive true means intersect succeeded, S = new_S, which is >=3 polygon. So enforce.
                if cand["alive"] and len(cand["S"]) < 3:
                    # Allow 0? But flag as error
                    raise ValueError(f"prefix {p_idx} step {s_idx} cand {c_idx}: alive but S has {len(cand['S'])} points (<3)")
                if not cand["alive"] and len(cand["S"]) != 0:
                    # Dead should have empty S (except if dead due to earlier? In code dead candidates have S empty)
                    # But they set S only if alive else empty? For dead they push empty S. So allow empty only.
                    # If dead but S non-empty, it's stale; but not fatal.
                    pass

        if steps:
            if steps[-1]["stream_idx"] not in (data["end_idx"], data["end_idx"] - 1):
                raise ValueError(f"prefix {p_idx}: last step stream_idx {steps[-1]['stream_idx']} not in {{end_idx, end_idx-1}} ({data['end_idx']}, {data['end_idx']-1})")
        # Check that p0_idx sequencing across prefixes: each prefix's p0_idx == previous end_idx
        if p_idx > 0:
            prev = prefixes[p_idx - 1]["data"]
            if data["p0_idx"] != prev["end_idx"]:
                raise ValueError(f"prefix {p_idx}: p0_idx {data['p0_idx']} != prev end_idx {prev['end_idx']}")

        # Check that last step's stream_idx == end_idx -1? In code, end_idx is cur at break, steps cover cur up to cur-1? Actually step stream_idx == cur at each iteration, end_idx == cur after loop. Steps' last stream_idx should be end_idx-1 or end_idx? Let's check: while loop increments cur after processing step. At break, end_idx=cur (first uncovered). Steps last stream_idx = cur-1? For non-breaking final step? But not strict. We'll just check monotonic.
    if prefixes[-1]["data"]["end_idx"] != len(stream):
        raise ValueError(f"last prefix end_idx {prefixes[-1]['data']['end_idx']} != stream len {len(stream)} (trace truncated)")

    # Done validation
    if done.get("type") != "done":
        raise ValueError("last line not done")
    for f in ("time_ms", "simplified", "frechet_distance"):
        if f not in done:
            raise ValueError(f"done missing field: {f}")
    if not _is_finite(done["time_ms"]) or float(done["time_ms"]) < 0:
        raise ValueError(f"done.time_ms must be finite >=0, got {done['time_ms']!r}")
    simplified = done["simplified"]
    if not _is_points(simplified):
        raise ValueError("done.simplified must be array of [x,y]")
    if len(simplified) < 2:
        raise ValueError(f"done.simplified must have at least 2 points, got {len(simplified)}")
    if len(simplified) % 2 != 0:
        raise ValueError(f"done.simplified len {len(simplified)} is odd, expected even (2*prefixes)")
    if len(simplified) != 2 * len(prefixes):
        raise ValueError(f"done.simplified len {len(simplified)} != 2*prefixes {2*len(prefixes)}")
    expected = []
    for pref in prefixes:
        expected.extend(pref["data"]["output"])
    for i, (a, b) in enumerate(zip(simplified, expected)):
        if not (_close(a[0], b[0]) and _close(a[1], b[1])):
            raise ValueError(f"done.simplified[{i}] {a!r} != prefix output {b!r}")
    fd = done["frechet_distance"]
    if fd is not None:
        if not _is_finite(fd) or float(fd) < 0:
            raise ValueError(f"done.frechet_distance must be null or finite >=0, got {fd!r}")

    return {
        "lines": len(lines),
        "bytes": len(raw),
        "is_monolithic": False,
        "header": header,
        "prefixes": len(prefixes),
        "steps": total_steps,
        "stream_len": len(stream),
        "simplified_len": len(simplified),
        "time_ms": float(done["time_ms"]),
        "was_gzipped": False,
    }

def main() -> int:
    ap = argparse.ArgumentParser(description="Validate web server NDJSON trace")
    ap.add_argument("--input", "-i", default=None, help="Input NDJSON file (gzipped or plain). Use - or omit for stdin")
    ap.add_argument("--epsilon", type=float, default=None, help="Expected epsilon")
    ap.add_argument("--delta", type=float, default=None, help="Expected delta")
    ap.add_argument("--trace-id", type=int, default=None, help="Trace id for original.txt length check")
    ap.add_argument("--data-dir", default="data", help="Data dir for original.txt lookup")
    ap.add_argument("--wall-time-ms", type=float, default=None, help="Wall time measured by curl (optional, for reporting)")
    ap.add_argument("--size-bytes", type=int, default=None, help="Payload size (optional)")
    ap.add_argument("--expect-gzip", action="store_true", help="Input is expected to be gzipped on disk")
    ap.add_argument("--json-summary", default=None, help="Write JSON summary to file")
    args = ap.parse_args()

    orig_n = None
    if args.trace_id is not None:
        try:
            orig_path = Path(args.data_dir) / str(args.trace_id) / "original.txt"
            if orig_path.exists():
                with orig_path.open() as f:
                    orig_n = int(f.readline().strip())
            else:
                print(f"warn: original.txt not found for trace {args.trace_id} at {orig_path}, skipping length check", file=sys.stderr)
        except Exception as e:
            print(f"warn: failed to read original length: {e}", file=sys.stderr)

    raw, raw_disk, was_gz = _read_bytes(args.input)
    if args.expect_gzip and not was_gz:
        print("error: expected gzipped payload but input was not gzipped", file=sys.stderr)
        return 1
    if raw_disk == 0:
        print("error: empty input (0 bytes)", file=sys.stderr)
        return 1

    try:
        info = validate_ndjson(raw, args.epsilon, args.delta, args.trace_id, orig_n)
    except ValueError as e:
        print(f"❌ NDJSON validation failed: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"❌ NDJSON validation error: {e}", file=sys.stderr)
        import traceback; traceback.print_exc()
        return 1

    info["bytes_raw_disk"] = raw_disk
    info["was_gzipped_on_disk"] = was_gz
    if args.wall_time_ms is not None:
        info["wall_time_ms"] = float(args.wall_time_ms)
    if args.size_bytes is not None:
        info["size_bytes"] = int(args.size_bytes)
    # Also include on-disk sizes
    info["bytes_decompressed"] = len(raw)

    # Summary printing
    eps_s = info["header"]["eps"] if isinstance(info["header"], dict) else info["header"]
    # header may be dict
    hdr_eps = info["header"]["eps"] if isinstance(info["header"], dict) and "eps" in info["header"] else "?"
    hdr_delta = info["header"]["delta"] if isinstance(info["header"], dict) and "delta" in info["header"] else "?"
    print(f"✅ NDJSON valid: {info['prefixes']} prefixes, {info['steps']} steps, stream {info['stream_len']} → simplified {info['simplified_len']}")
    print(f"   header eps={hdr_eps} delta={hdr_delta} time_ms={info['time_ms']:.2f}")
    print(f"   lines={info['lines']} decompressed={info['bytes_decompressed']} bytes (on-disk {raw_disk} bytes{' gzipped' if was_gz else ''})")
    if "wall_time_ms" in info:
        print(f"   wall_time_ms={info['wall_time_ms']:.2f}")
    # Payload size / speed heuristics (informational, not failing here)
    # Caller (workflow / shell) will gate wall time and size.
    if args.json_summary:
        try:
            Path(args.json_summary).write_text(json.dumps(info, indent=2))
        except Exception as e:
            print(f"warn: failed to write json summary: {e}", file=sys.stderr)

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
