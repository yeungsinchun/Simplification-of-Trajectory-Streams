#!/usr/bin/env python3
"""
Flask backend for the trajectory simplification web visualizer.

Accepts an uploaded original.txt trajectory file + epsilon/delta params,
runs the C++ simplify binary with --web-server --json-stream, and streams
the trace as NDJSON (header, one prefix per line, done).
"""
import os
import subprocess
import time
import uuid
import shutil
import json
import gzip
import re
from pathlib import Path
from flask import Flask, request, jsonify, send_from_directory, Response

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024  # 50MB max upload

# Resolve repo root (server.py lives in web/, repo root is parent)
REPO_ROOT = Path(__file__).parent.parent.resolve()
SIMPLIFY_BIN = REPO_ROOT / "build" / "simplify"
DOTS_BIN = REPO_ROOT / "build" / "dots"
DP_BIN = REPO_ROOT / "build" / "dp"
SQUISH_BIN = REPO_ROOT / "build" / "squish"
FRECHET_BIN  = REPO_ROOT / "scripts" / "frechet.jl"
DATA_DIR = REPO_ROOT / "data"
WEB_DIR = REPO_ROOT / "web"

# Baseline polyline files under data/<id>/ (first hit wins per layer key).
LAYER_FILES = {
    'dots': ('dots_simplified.txt', 'DOTS.txt'),
    'dp': ('dp_simplified.txt', 'DP.txt'),
    'squish': ('squish_simplified.txt', 'SQUISH.txt'),
    'simplify': ('simplify.txt',),
}

BASELINE_ALGOS = {
    'dots': {
        'label': 'DOTS',
        'bin': DOTS_BIN,
        'core_ms_re': r'DOTS_CORE_MS:\s*([\d.]+)',
        'files': LAYER_FILES['dots'],
    },
    'dp': {
        'label': 'DP',
        'bin': DP_BIN,
        'core_ms_re': r'DP_CORE_MS:\s*([\d.]+)',
        'files': LAYER_FILES['dp'],
    },
    'squish': {
        'label': 'SQUISH',
        'bin': SQUISH_BIN,
        'core_ms_re': r'SQUISH_CORE_MS:\s*([\d.]+)',
        'files': LAYER_FILES['squish'],
    },
}


def read_xy_polyline(path: Path):
    """Parse N\\n x y files used by simplify/dots. Returns list[[x,y]] or None."""
    if not path.exists():
        return None
    try:
        with open(path, 'r') as f:
            first = f.readline().strip()
            if not first:
                return None
            n = int(first)
            points = []
            for _ in range(n):
                line = f.readline()
                if not line:
                    break
                parts = line.split()
                if len(parts) >= 2:
                    points.append([float(parts[0]), float(parts[1])])
            return points if points else None
    except Exception:
        return None


def load_trace_layers(trace_id: int):
    """Load available compare layers for a preloaded trace id."""
    trace_dir = DATA_DIR / str(trace_id)
    layers = {}
    sources = {}
    for key, names in LAYER_FILES.items():
        for name in names:
            pts = read_xy_polyline(trace_dir / name)
            if pts is not None:
                layers[key] = pts
                sources[key] = name
                break
    original = read_xy_polyline(trace_dir / 'original.txt')
    if original is not None:
        layers['original'] = original
        sources['original'] = 'original.txt'
    return layers, sources


def trace_json_response(tracejson):
    """Return compact JSON, gzipped when the client supports it."""
    payload = json.dumps(tracejson, separators=(',', ':')).encode('utf-8')
    response = Response(payload, mimetype='application/json')

    # Large traces can exceed Cloud Run's 32 MiB HTTP/1 response limit.
    # Browsers advertise gzip support, and Cloud Run counts the compressed
    # bytes sent by the container.
    accepted_encodings = request.headers.get('Accept-Encoding', '').lower()
    if 'gzip' in accepted_encodings:
        response.set_data(gzip.compress(payload, compresslevel=6))
        response.headers['Content-Encoding'] = 'gzip'
        response.headers['Vary'] = 'Accept-Encoding'

    return response


def stream_simplify_trace(cmd, label):
    """Run simplify with --json-stream and yield NDJSON lines as they are produced."""
    def generate():
        # Keep stderr on a pipe but drain it in a thread so a full stderr
        # buffer cannot deadlock the C++ process while we read stdout.
        proc = subprocess.Popen(
            cmd,
            cwd=REPO_ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )
        stderr_chunks = []

        def _drain_stderr():
            assert proc.stderr is not None
            for chunk in proc.stderr:
                stderr_chunks.append(chunk)

        import threading
        stderr_thread = threading.Thread(target=_drain_stderr, daemon=True)
        stderr_thread.start()

        emitted = False
        try:
            assert proc.stdout is not None
            for line in proc.stdout:
                if line:
                    emitted = True
                    yield line if line.endswith('\n') else line + '\n'
            proc.wait(timeout=300)
            stderr_thread.join(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
            yield json.dumps({'type': 'error', 'message': 'Processing timeout (>5min)'}) + '\n'
            return
        except Exception as e:
            proc.kill()
            proc.wait()
            yield json.dumps({'type': 'error', 'message': str(e)}) + '\n'
            return

        err = ''.join(stderr_chunks).strip()
        if err:
            print(f"[{label}] C++ stderr:\n{err}")
        if proc.returncode != 0:
            message = err or 'Binary execution failed'
            print(f"[{label}] simplify failed: {message}")
            yield json.dumps({'type': 'error', 'message': message}) + '\n'

    return Response(
        generate(),
        mimetype='application/x-ndjson',
        headers={
            'Cache-Control': 'no-cache',
            'X-Accel-Buffering': 'no',
            # Disable Flask/Werkzeug buffering hints for proxies.
            'Connection': 'keep-alive',
        },
    )


# --- Static file serving ---

@app.route('/')
def index():
    return send_from_directory(WEB_DIR, 'index.html')

@app.route('/<path:filename>')
def static_files(filename):
    return send_from_directory(WEB_DIR, filename)

# --- Preloaded traces endpoint ---

@app.route('/api/traces', methods=['GET'])
def list_traces():
    """
    GET /api/traces
    Returns a list of available preloaded trace IDs.
    Sample traces (51-53) are returned first with their label and recommended
    epsilon/delta pulled from meta.json; regular traces (1-50) follow.
    """
    def read_trace_entry(i):
        trace_dir = DATA_DIR / str(i)
        orig = trace_dir / "original.txt"
        if not (trace_dir.exists() and orig.exists()):
            return None
        entry = {'id': i, 'n_points': None, 'label': None,
                 'epsilon': None, 'delta': None}
        try:
            with open(orig, 'r') as f:
                entry['n_points'] = int(f.readline().strip())
        except Exception:
            pass
        meta_path = trace_dir / "meta.json"
        if meta_path.exists():
            try:
                with open(meta_path, 'r') as f:
                    meta = json.load(f)
                entry['label']   = meta.get('label')
                entry['epsilon'] = meta.get('epsilon')
                entry['delta']   = meta.get('delta')
            except Exception:
                pass
        return entry

    traces = []
    # Sample traces first (51-53)
    for i in (51, 52, 53):
        e = read_trace_entry(i)
        if e:
            traces.append(e)
    # Regular traces
    for i in range(1, 51):
        e = read_trace_entry(i)
        if e:
            traces.append(e)
    return jsonify({'traces': traces})

@app.route('/api/trace/<int:trace_id>/original', methods=['GET'])
def get_original_trace(trace_id):
    """
    GET /api/trace/<id>/original
    Returns the original trajectory points only (for quick preview).
    """
    trace_dir = DATA_DIR / str(trace_id)
    original_file = trace_dir / "original.txt"
    if not original_file.exists():
        return jsonify({'error': f'Trace {trace_id} not found'}), 404
    
    try:
        with open(original_file, 'r') as f:
            lines = f.readlines()
            n = int(lines[0].strip())
            points = []
            for i in range(1, min(n + 1, len(lines))):
                parts = lines[i].strip().split()
                if len(parts) >= 2:
                    points.append([float(parts[0]), float(parts[1])])
        
        return jsonify({'points': points, 'trace_id': trace_id})
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/trace/<int:trace_id>/compare', methods=['GET'])
def get_trace_compare(trace_id):
    """
    GET /api/trace/<id>/compare
    Query:
      algorithm: none|dots|dp|squish (default none)
      run: 1 to execute the baseline binary
      lssd: DOTS LSSD threshold (default 1e6 / 1000K)
      epsilon: DP PED epsilon (default 0.9)
      ratio: SQUISH compression ratio in (0, 1] (default 0.2)
    """
    trace_dir = DATA_DIR / str(trace_id)
    if not trace_dir.exists() or not (trace_dir / 'original.txt').exists():
        return jsonify({'error': f'Trace {trace_id} not found'}), 404

    algorithm = (request.args.get('algorithm') or 'none').lower().strip()
    if algorithm in ('', 'none', 'null'):
        algorithm = 'none'
    run = request.args.get('run', '0') in ('1', 'true', 'yes')

    if algorithm not in ('none',) and algorithm not in BASELINE_ALGOS:
        return jsonify({'error': f'Unknown algorithm: {algorithm}'}), 400

    lssd = 1e6
    dp_eps = 0.9
    squish_ratio = 0.2
    try:
        if request.args.get('lssd') is not None:
            lssd = float(request.args.get('lssd'))
            if lssd <= 0:
                return jsonify({'error': 'lssd must be positive'}), 400
        if request.args.get('epsilon') is not None:
            dp_eps = float(request.args.get('epsilon'))
            if dp_eps <= 0:
                return jsonify({'error': 'epsilon must be positive'}), 400
        if request.args.get('ratio') is not None:
            squish_ratio = float(request.args.get('ratio'))
            if squish_ratio <= 0 or squish_ratio > 1:
                return jsonify({'error': 'ratio must be in (0, 1]'}), 400
    except ValueError:
        return jsonify({'error': 'Invalid numeric baseline parameter'}), 400

    baseline_core_ms = None
    baseline_error = None
    baseline_label = BASELINE_ALGOS[algorithm]['label'] if algorithm in BASELINE_ALGOS else None

    if algorithm in BASELINE_ALGOS and run:
        meta = BASELINE_ALGOS[algorithm]
        binary = meta['bin']
        if not binary.exists():
            baseline_error = f'{algorithm} binary not found at {binary}'
        else:
            original = trace_dir / 'original.txt'
            out_name = meta['files'][0]
            out_path = trace_dir / out_name
            if algorithm == 'dots':
                cmd = [str(binary), str(trace_id), '-lssd', str(lssd)]
                cwd = REPO_ROOT / 'build'
            elif algorithm == 'dp':
                cmd = [str(binary), str(original), str(dp_eps), str(out_path)]
                cwd = REPO_ROOT
            else:
                cmd = [str(binary), str(original), str(squish_ratio), str(out_path)]
                cwd = REPO_ROOT
            try:
                started = time.perf_counter()
                result = subprocess.run(
                    cmd,
                    cwd=cwd,
                    capture_output=True,
                    text=True,
                    timeout=120,
                )
                elapsed_ms = (time.perf_counter() - started) * 1000.0
                out = (result.stdout or '') + (result.stderr or '')
                if result.returncode != 0:
                    baseline_error = (result.stderr or result.stdout or f'{algorithm} failed').strip()
                else:
                    m = re.search(meta['core_ms_re'], out)
                    baseline_core_ms = float(m.group(1)) if m else elapsed_ms
            except Exception as e:
                baseline_error = str(e)

    layers, sources = load_trace_layers(trace_id)
    if algorithm in BASELINE_ALGOS and run and baseline_error:
        layers.pop(algorithm, None)
        sources.pop(algorithm, None)
    elif algorithm in BASELINE_ALGOS and algorithm in layers:
        layers['baseline'] = layers[algorithm]
        sources['baseline'] = sources.get(algorithm)

    metrics = {
        'algorithm': algorithm,
        'baseline_label': baseline_label,
        'lssd': lssd if algorithm == 'dots' else None,
        'epsilon': dp_eps if algorithm == 'dp' else None,
        'ratio': squish_ratio if algorithm == 'squish' else None,
        'original_points': len(layers['original']) if 'original' in layers else None,
        'simplify_points': len(layers['simplify']) if 'simplify' in layers else None,
        'baseline_points': len(layers['baseline']) if 'baseline' in layers else None,
        'baseline_core_ms': baseline_core_ms,
        'simplify_core_ms': None,
        'dots_points': len(layers['dots']) if 'dots' in layers else None,
        'dots_core_ms': baseline_core_ms if algorithm == 'dots' else None,
    }
    return jsonify({
        'trace_id': trace_id,
        'algorithm': algorithm,
        'layers': layers,
        'sources': sources,
        'available': sorted(layers.keys()),
        'metrics': metrics,
        'baseline_error': baseline_error,
        'dots_error': baseline_error if algorithm == 'dots' else None,
        'algorithms': [
            {'id': 'dots', 'label': 'DOTS', 'params': [
                {'id': 'lssd', 'label': 'DOTS LSSD', 'type': 'number', 'default': 1e6, 'unit': 'K', 'min': 1e-3, 'step': 1},
            ]},
            {'id': 'dp', 'label': 'DP', 'params': [
                {'id': 'epsilon', 'label': 'DP PED ε', 'type': 'number', 'default': 0.9, 'min': 1e-9, 'step': 0.1},
            ]},
            {'id': 'squish', 'label': 'SQUISH', 'params': [
                {'id': 'ratio', 'label': 'SQUISH Ratio', 'type': 'number', 'default': 20, 'unit': '%', 'min': 0.01, 'max': 100, 'step': 1},
            ]},
        ],
    })


@app.route('/api/trace/<int:trace_id>/frechet/<curve>', methods=['GET'])
def frechet_existing_curve(trace_id, curve):
    """
    GET /api/trace/<id>/frechet/<curve>
    Frechet distance between original.txt and an existing result file.
    curve: simplify | dots
    """
    curve = curve.lower()
    name_map = {
        'simplify': ('simplify.txt',),
        'dots': ('dots_simplified.txt', 'DOTS.txt'),
        'dp': ('dp_simplified.txt', 'DP.txt'),
        'squish': ('squish_simplified.txt', 'SQUISH.txt'),
    }
    if curve not in name_map:
        return jsonify({'error': f'Unknown curve {curve}'}), 400

    trace_dir = DATA_DIR / str(trace_id)
    original = trace_dir / 'original.txt'
    if not original.exists():
        return jsonify({'error': f'Trace {trace_id} not found'}), 404

    simplified = None
    for name in name_map[curve]:
        candidate = trace_dir / name
        if candidate.exists():
            simplified = candidate
            break
    if simplified is None:
        return jsonify({'error': f'No {curve} polyline file for trace {trace_id}'}), 404

    request_id = f'{trace_id}_{curve}_{uuid.uuid4().hex[:8]}'
    try:
        distance = _run_frechet(request_id, original, simplified)
        return jsonify({
            'trace_id': trace_id,
            'curve': curve,
            'path': simplified.name,
            'distance': distance,
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/trace/<int:trace_id>', methods=['GET'])
def get_trace(trace_id):
    """
    GET /api/trace/<id>
    Query params:
      - epsilon: float (default 0.5)
      - delta: float (default 300)
    
    Returns JSON trace for preloaded trace ID.
    """
    trace_dir = DATA_DIR / str(trace_id)
    if not trace_dir.exists() or not (trace_dir / "original.txt").exists():
        return jsonify({'error': f'Trace {trace_id} not found'}), 404
    
    try:
        epsilon = float(request.args.get('epsilon', 0.5))
        delta = float(request.args.get('delta', 300))
    except ValueError:
        return jsonify({'error': 'Invalid epsilon or delta'}), 400
    
    try:
        import time

        cmd = [
            str(SIMPLIFY_BIN),
            '--in', str(trace_id),
            '-e', str(epsilon),
            '-d', str(delta),
            '--web-server',
            '--json-stream',
        ]

        print(f"[Trace {trace_id}] Streaming simplify with eps={epsilon}, delta={delta}")
        start_time = time.time()
        response = stream_simplify_trace(cmd, f"Trace {trace_id}")

        @response.call_on_close
        def _log_done():
            elapsed = time.time() - start_time
            print(f"[Trace {trace_id}] Stream finished in {elapsed:.2f}s")

        return response
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# --- Trace generation endpoint ---

@app.route('/api/trace', methods=['POST'])
def generate_trace():
    """
    POST /api/trace
    Form data:
      - file: the trajectory file (N\\nx1 y1\\n... format)
      - epsilon: float (default 0.5)
      - delta: float (default 300)
    
    Returns JSON trace or error.
    """
    if 'file' not in request.files:
        return jsonify({'error': 'No file provided'}), 400
    
    uploaded_file = request.files['file']
    if uploaded_file.filename == '':
        return jsonify({'error': 'Empty filename'}), 400
    
    try:
        epsilon = float(request.form.get('epsilon', 0.5))
        delta = float(request.form.get('delta', 300))
    except ValueError:
        return jsonify({'error': 'Invalid epsilon or delta'}), 400
    
    # Validate basic format (first line should be an integer N)
    content = uploaded_file.read().decode('utf-8')
    lines = content.strip().split('\n')
    if len(lines) < 1:
        return jsonify({'error': 'Empty file'}), 400
    try:
        n = int(lines[0])
        if n < 2 or len(lines) < n + 1:
            return jsonify({'error': f'Invalid format: expected {n} points'}), 400
    except ValueError:
        return jsonify({'error': 'First line must be point count'}), 400
    
    # Create a temp directory for this request. The simplify binary's --in
    # flag parses its argument with std::stoi, so the id must be numeric
    # (it builds the path data/<id>/original.txt internally).
    temp_id = str(900000000 + uuid.uuid4().int % 99999999)
    temp_dir = DATA_DIR / temp_id
    temp_dir.mkdir(parents=True, exist_ok=True)
    temp_file = temp_dir / "original.txt"
    
    try:
        # Write uploaded content
        with open(temp_file, 'w') as f:
            f.write(content)
        
        import time

        cmd = [
            str(SIMPLIFY_BIN),
            '--in', temp_id,
            '-e', str(epsilon),
            '-d', str(delta),
            '--web-server',
            '--json-stream',
        ]

        print(f"[Upload] Streaming simplify with eps={epsilon}, delta={delta}")
        start_time = time.time()
        response = stream_simplify_trace(cmd, "Upload")

        @response.call_on_close
        def _cleanup_upload():
            elapsed = time.time() - start_time
            print(f"[Upload] Stream finished in {elapsed:.2f}s")
            if temp_dir.exists():
                shutil.rmtree(temp_dir, ignore_errors=True)

        return response
    
    except Exception as e:
        if temp_dir.exists():
            shutil.rmtree(temp_dir, ignore_errors=True)
        return jsonify({'error': str(e)}), 500


# --- Fréchet distance ---

def _run_frechet(request_id: str, original_path: Path, simplified_path: Path):
    """Run Julia and return the Fréchet distance."""
    print(f"[Fréchet {request_id}] Starting computation...")
    tmp_id = str(800000000 + uuid.uuid4().int % 99999999)
    tmp_dir = DATA_DIR / tmp_id
    try:
        # Use batch mode with a non-zero ID to properly trigger batch mode
        tmp_dir.mkdir(parents=True, exist_ok=True)
        tmp_orig = tmp_dir / "original.txt"
        shutil.copy2(original_path, tmp_orig)

        cmd = ["julia", str(FRECHET_BIN), "--raw",
               "--id", str(tmp_id),
               "--batch", str(simplified_path)]

        print(f"[Fréchet {request_id}] Running: {' '.join(cmd)}")
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=120
        )

        if result.returncode != 0:
            print(f"[Fréchet {request_id}] Failed with exit code {result.returncode}")
            print(f"[Fréchet {request_id}] stderr: {result.stderr}")
            raise RuntimeError(result.stderr.strip() or f"exit {result.returncode}")

        # Output is "<basename>: <distance>" in batch mode
        out = result.stdout.strip()
        print(f"[Fréchet {request_id}] Output: {out}")
        # Parse: "simplify.txt: 123.45" -> 123.45
        if ":" in out:
            dist_str = out.split(":")[-1].strip()
        else:
            dist_str = out.split()[-1].rstrip(":,")
        distance = float(dist_str)

        print(f"[Fréchet {request_id}] Success! Distance = {distance}")
        return distance
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


@app.route('/api/frechet', methods=['POST'])
def start_frechet():
    """
    POST /api/frechet
    JSON body:
      { "trace_id": <int>, "eps": <float>, "delta": <float> }
    or
      { "file_content": "<str>", "eps": <float>, "delta": <float> }

    Computes and returns { "distance": <float> } in the same request so the
    result cannot be lost when Cloud Run routes requests across instances.
    """
    data = request.get_json(force=True, silent=True) or {}
    eps   = float(data.get('eps', 0.5))
    delta = float(data.get('delta', 300))

    request_id = str(uuid.uuid4())
    print(f"[Fréchet] POST /api/frechet - request_id={request_id}, eps={eps}, delta={delta}")

    # --- Obtain original trajectory and build simplified curve ---
    tmp_root = DATA_DIR / f"frechet_{request_id}"
    tmp_root.mkdir(parents=True, exist_ok=True)
    original_path   = tmp_root / "original.txt"
    simplified_path = tmp_root / "simplify.txt"
    tmp_simp_dir = None

    try:
        if 'trace_id' in data:
            trace_dir = DATA_DIR / str(int(data['trace_id']))
            src_orig = trace_dir / "original.txt"
            print(f"[Fréchet {request_id}] Using trace_id={data['trace_id']}, path={src_orig}")
            if not src_orig.exists():
                return jsonify({'error': f"Trace {data['trace_id']} not found"}), 404
            shutil.copy2(src_orig, original_path)
        elif 'file_content' in data:
            print(f"[Fréchet {request_id}] Using file_content")
            original_path.write_text(data['file_content'])
        else:
            return jsonify({'error': 'Provide trace_id or file_content'}), 400

        # Run simplify to get the simplified curve
        tmp_id = str(900000000 + uuid.uuid4().int % 99999999)
        tmp_simp_dir = DATA_DIR / tmp_id
        tmp_simp_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(original_path, tmp_simp_dir / "original.txt")
        cmd = [str(SIMPLIFY_BIN), '--in', tmp_id, '--out', '-e', str(eps), '-d', str(delta)]
        print(f"[Fréchet {request_id}] Running simplify: {' '.join(cmd)}")
        result = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True, timeout=300)
        # simplify writes data/<id>/simplify.txt
        simp_out = tmp_simp_dir / "simplify.txt"
        if result.returncode != 0 or not simp_out.exists():
            print(f"[Fréchet {request_id}] Simplify failed: {result.stderr}")
            return jsonify({'error': 'Simplify failed', 'stderr': result.stderr}), 500
        shutil.copy2(simp_out, simplified_path)
        shutil.rmtree(tmp_simp_dir, ignore_errors=True)
        tmp_simp_dir = None
        print(f"[Fréchet {request_id}] Simplify succeeded")

        distance = _run_frechet(request_id, original_path, simplified_path)
        return jsonify({'distance': distance})
    except Exception as e:
        print(f"[Fréchet {request_id}] Exception: {e}")
        return jsonify({'error': str(e)}), 500
    finally:
        if tmp_simp_dir is not None:
            shutil.rmtree(tmp_simp_dir, ignore_errors=True)
        shutil.rmtree(tmp_root, ignore_errors=True)


if __name__ == '__main__':
    # Check that binary exists
    if not SIMPLIFY_BIN.exists():
        print(f"Error: simplify binary not found at {SIMPLIFY_BIN}")
        print("Build it first: cmake -B build && cmake --build build")
        exit(1)
    
    port = int(os.environ.get('PORT', 5050))
    # Bind to 0.0.0.0 by default. Cloud Run's request router and Docker's
    # published-port forwarding both require a non-loopback bind inside the
    # container; binding to 127.0.0.1 works for `python server.py` on a dev
    # machine but breaks in any containerized deployment. Set HOST=127.0.0.1
    # explicitly if you need loopback-only on a dev box.
    host = os.environ.get('HOST', '0.0.0.0')
    print(f"Repo root: {REPO_ROOT}")
    print(f"Binary: {SIMPLIFY_BIN}")
    print(f"Starting Flask server on http://{host}:{port}")
    # Note: on macOS, port 5000 is claimed by AirPlay Receiver (AirTunes),
    # which silently returns 403 for any request. Use 5050 or set PORT= to
    # override.
    # threaded=True so concurrent requests (e.g. the long-running Julia
    # frechet job) don't block simpler endpoints.
    app.run(debug=False, host=host, port=port, threaded=True)
