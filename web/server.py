#!/usr/bin/env python3
"""
Flask backend for the trajectory simplification web visualizer.

Accepts an uploaded original.txt trajectory file + epsilon/delta params,
runs the C++ simplify binary with --web-server --json-stream, and streams
the trace as NDJSON (header, one prefix per line, done).
"""
import os
import subprocess
import uuid
import shutil
import json
import gzip
from pathlib import Path
from flask import Flask, request, jsonify, send_from_directory, Response

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024  # 50MB max upload

# Resolve repo root (server.py lives in web/, repo root is parent)
REPO_ROOT = Path(__file__).parent.parent.resolve()
SIMPLIFY_BIN = REPO_ROOT / "build" / "simplify"
FRECHET_BIN  = REPO_ROOT / "scripts" / "frechet.jl"
DATA_DIR = REPO_ROOT / "data"
WEB_DIR = REPO_ROOT / "web"

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
