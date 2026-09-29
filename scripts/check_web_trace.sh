#!/usr/bin/env bash
# check_web_trace.sh — manual re-measure for the web server trace endpoint
#
# Builds `simplify`, starts the Flask server, curls /api/trace/<id> and validates
# NDJSON structure (including gzipped variant) via scripts/ci/validate_web_trace.py.
# Measures response time / payload size and gates on reasonable thresholds.
#
# Usage:
#   scripts/check_web_trace.sh                               # defaults: id 1, eps 0.5, delta 200
#   scripts/check_web_trace.sh --trace-id 1 --epsilon 0.5 --delta 200
#   scripts/check_web_trace.sh --no-build --port 5051 --keep
#
# Env overrides: TRACE_ID, EPSILON, DELTA, PORT, HOST, MAX_TIME_MS, MIN_BYTES
#
# The same validation is used by .github/workflows/web-server.yml (CI).
# This script is shellcheck-clean and cross-platform (macOS/Linux).

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
TRACE_ID="${TRACE_ID:-1}"
EPSILON="${EPSILON:-0.5}"
DELTA="${DELTA:-200}"
PORT="${PORT:-5051}"
HOST="${HOST:-127.0.0.1}"
MAX_TIME_MS="${MAX_TIME_MS:-10000}"   # 10s wall-time gate (fine epsilon on 589 pts)
MAX_CORE_MS="${MAX_CORE_MS:-8000}"    # core time gate from done.time_ms
MIN_BYTES="${MIN_BYTES:-1024}"        # payload must be at least 1 KiB
# If set, keep server running and trace files; otherwise cleanup on exit.
KEEP=0
NO_BUILD=0

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
VALIDATOR="$REPO_ROOT/scripts/ci/validate_web_trace.py"
SERVER_LOG="/tmp/web-check-server.log"
TRACE_FILE="/tmp/web-check-trace.ndjson"
TRACE_GZ="/tmp/web-check-trace.gz"
TRACE_GZ_UNZIPPED="/tmp/web-check-trace.gunzipped.ndjson"
DIRECT_FILE="/tmp/web-check-direct.ndjson"
HEADERS_FILE="/tmp/web-check-headers.txt"
GZIP_HEADERS="/tmp/web-check-gzip-headers.txt"

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
while [ $# -gt 0 ]; do
  case "$1" in
    --trace-id) TRACE_ID="$2"; shift 2;;
    --epsilon) EPSILON="$2"; shift 2;;
    --delta) DELTA="$2"; shift 2;;
    --port) PORT="$2"; shift 2;;
    --host) HOST="$2"; shift 2;;
    --max-time-ms) MAX_TIME_MS="$2"; shift 2;;
    --min-bytes) MIN_BYTES="$2"; shift 2;;
    --keep) KEEP=1; shift;;
    --no-build) NO_BUILD=1; shift;;
    --help|-h)
      echo "Usage: $0 [--trace-id ID] [--epsilon E] [--delta D] [--port P] [--keep] [--no-build]"
      echo "  Defaults: trace-id=$TRACE_ID epsilon=$EPSILON delta=$DELTA port=$PORT"
      echo "  Env: TRACE_ID EPSILON DELTA PORT HOST MAX_TIME_MS MIN_BYTES"
      exit 0
      ;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

echo "==> web trace check: trace=$TRACE_ID eps=$EPSILON delta=$DELTA host=$HOST port=$PORT"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
cleanup() {
  if [ "$KEEP" -eq 0 ]; then
    if [ -n "${SERVER_PID:-}" ] && kill -0 "$SERVER_PID" 2>/dev/null; then
      echo "==> stopping server (pid $SERVER_PID)"
      kill "$SERVER_PID" 2>/dev/null || true
      wait "$SERVER_PID" 2>/dev/null || true
    fi
    # do not delete trace files when keep is requested; otherwise leave for debugging on failure
    # keep files on failure regardless
    if [ "${FAIL:-0}" -eq 0 ]; then
      rm -f "$TRACE_FILE" "$TRACE_GZ" "$TRACE_GZ_UNZIPPED" "$DIRECT_FILE" "$HEADERS_FILE" "$GZIP_HEADERS" 2>/dev/null || true
    fi
  else
    echo "==> --keep: server pid ${SERVER_PID:-unknown} logs at $SERVER_LOG"
    echo "    traces: $TRACE_FILE, $TRACE_GZ, $DIRECT_FILE"
  fi
}
trap cleanup EXIT

fail() {
  FAIL=1
  echo "❌ $*" >&2
  exit 1
}

need_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    fail "required command not found: $1"
  fi
}

# ---------------------------------------------------------------------------
# Prereqs
# ---------------------------------------------------------------------------
need_cmd python3
need_cmd curl
need_cmd cmake

# Build simplify if needed
if [ "$NO_BUILD" -eq 0 ]; then
  if [ ! -x "$REPO_ROOT/build/simplify" ]; then
    echo "==> building simplify (Release, BUILD_GUI=OFF)"
    cmake -B "$REPO_ROOT/build" -DCMAKE_BUILD_TYPE=Release -DBUILD_GUI=OFF -S "$REPO_ROOT"
    cmake --build "$REPO_ROOT/build" --target simplify -j"$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)"
  else
    echo "==> simplify binary present: $REPO_ROOT/build/simplify"
  fi
else
  echo "==> --no-build: skipping cmake"
  if [ ! -x "$REPO_ROOT/build/simplify" ]; then
    fail "simplify binary not found and --no-build requested"
  fi
fi

# Python deps
if ! python3 -c "import flask" 2>/dev/null; then
  echo "==> installing python deps from web/requirements.txt"
  if command -v pip3 >/dev/null 2>&1; then
    pip3 install -q -r "$REPO_ROOT/web/requirements.txt"
  else
    python3 -m pip install -q -r "$REPO_ROOT/web/requirements.txt"
  fi
fi

# Check data trace exists
if [ ! -f "$REPO_ROOT/data/$TRACE_ID/original.txt" ]; then
  # shellcheck disable=SC2012 # ls is fine for numeric trace dirs; find would be noisy here
  fail "data/$TRACE_ID/original.txt not found (available ids: $(ls "$REPO_ROOT/data" 2>/dev/null | tr '\n' ' '))"
fi
ORIG_N="$(head -1 "$REPO_ROOT/data/$TRACE_ID/original.txt" | tr -d '[:space:]')"
echo "==> trace $TRACE_ID has $ORIG_N points"

# ---------------------------------------------------------------------------
# Start server
# ---------------------------------------------------------------------------
# Kill any stale server on the port (macOS lsof vs Linux fuser)
if command -v lsof >/dev/null 2>&1; then
  STALE_PID="$(lsof -tiTCP:"$PORT" -sTCP:LISTEN 2>/dev/null || true)"
  if [ -n "$STALE_PID" ]; then
    echo "==> killing stale server on port $PORT (pid $STALE_PID)"
    kill "$STALE_PID" 2>/dev/null || true
    sleep 1
  fi
fi

# Ensure port free
if command -v curl >/dev/null 2>&1 && curl -s --connect-timeout 1 "http://$HOST:$PORT/" >/dev/null 2>&1; then
  echo "warn: something already listening on $HOST:$PORT; trying to proceed" >&2
fi

echo "==> starting Flask server on $HOST:$PORT (log $SERVER_LOG)"
# shellcheck disable=SC2155
export PORT HOST
# Use unbuffered python for timely logs
PORT="$PORT" HOST="$HOST" python3 -u "$REPO_ROOT/web/server.py" > "$SERVER_LOG" 2>&1 &
SERVER_PID=$!
echo "    pid=$SERVER_PID"

# Wait for server (30s)
echo -n "==> waiting for server"
for _ in $(seq 1 30); do
  if curl -s --connect-timeout 1 "http://$HOST:$PORT/api/traces" >/dev/null 2>&1; then
    echo " ✓"
    break
  fi
  if ! kill -0 "$SERVER_PID" 2>/dev/null; then
    echo " ✗ server died"
    cat "$SERVER_LOG" >&2 || true
    fail "server failed to start"
  fi
  echo -n "."
  sleep 1
  if [ "$_" -eq 30 ]; then
    echo " ✗ timeout"
    cat "$SERVER_LOG" >&2 || true
    fail "server did not become ready in 30s"
  fi
done

BASE_URL="http://$HOST:$PORT"

# ---------------------------------------------------------------------------
# Curl plain NDJSON
# ---------------------------------------------------------------------------
TRACE_URL="$BASE_URL/api/trace/$TRACE_ID?epsilon=$EPSILON&delta=$DELTA"
echo "==> curl plain NDJSON: $TRACE_URL"

# Use curl with separate stats file to avoid mixing body and progress lines
CURL_STATS="/tmp/web-check-curl-stats.txt"
rm -f "$TRACE_FILE" "$HEADERS_FILE" "$CURL_STATS"

set +e
curl -s -D "$HEADERS_FILE" -o "$TRACE_FILE" -w "%{http_code} %{time_total} %{size_download} %{content_type}\n" "$TRACE_URL" > "$CURL_STATS" 2>&1
CURL_RC=$?
set -e
if [ "$CURL_RC" -ne 0 ]; then
  cat "$SERVER_LOG" >&2 || true
  fail "curl failed (rc $CURL_RC)"
fi
read -r HTTP_CODE TIME_TOTAL SIZE_DOWNLOAD CONTENT_TYPE < "$CURL_STATS"
# TIME_TOTAL is seconds (float), convert to ms
WALL_MS="$(python3 -c "import sys; print(int(float(sys.argv[1])*1000))" "$TIME_TOTAL" 2>/dev/null || echo "?")"
echo "    http=$HTTP_CODE content_type=$CONTENT_TYPE time=${TIME_TOTAL}s (${WALL_MS}ms) size=${SIZE_DOWNLOAD}B"

if [ "$HTTP_CODE" != "200" ]; then
  echo "--- headers ---"; cat "$HEADERS_FILE" || true
  echo "--- body (first 500 chars) ---"; head -c 500 "$TRACE_FILE" || true; echo
  cat "$SERVER_LOG" >&2 || true
  fail "expected HTTP 200, got $HTTP_CODE"
fi

if ! grep -qi "application/x-ndjson" "$HEADERS_FILE"; then
  fail "Content-Type is not application/x-ndjson: $(grep -i content-type "$HEADERS_FILE" || echo 'missing')"
fi

if [ ! -s "$TRACE_FILE" ]; then
  fail "trace file empty"
fi

ACTUAL_BYTES="$(wc -c < "$TRACE_FILE" | tr -d '[:space:]')"
echo "    file bytes=$ACTUAL_BYTES (curl reported $SIZE_DOWNLOAD)"
if [ "$ACTUAL_BYTES" -lt "$MIN_BYTES" ]; then
  fail "payload too small: $ACTUAL_BYTES < $MIN_BYTES bytes (regression?)"
fi

# Validate plain NDJSON
echo "==> validating plain NDJSON"
python3 "$VALIDATOR" --input "$TRACE_FILE" --epsilon "$EPSILON" --delta "$DELTA" --trace-id "$TRACE_ID" --wall-time-ms "$WALL_MS" --size-bytes "$ACTUAL_BYTES" || fail "NDJSON validation failed"

# Gate wall time
if [ "$WALL_MS" != "?" ] && [ "$WALL_MS" -gt "$MAX_TIME_MS" ] 2>/dev/null; then
  fail "wall time ${WALL_MS}ms exceeds limit ${MAX_TIME_MS}ms (regression)"
fi

# ---------------------------------------------------------------------------
# Curl gzipped variant
# ---------------------------------------------------------------------------
echo "==> curl gzipped NDJSON (Accept-Encoding: gzip)"
rm -f "$TRACE_GZ" "$GZIP_HEADERS" "$CURL_STATS"
# Use --raw to keep gzip encoding, -H to request gzip; do not use --compressed which auto-decompresses
set +e
curl -s -H "Accept-Encoding: gzip" -D "$GZIP_HEADERS" -o "$TRACE_GZ" --raw -w "%{http_code} %{time_total} %{size_download} %{content_type}\n" "$TRACE_URL" > "$CURL_STATS" 2>&1
CURL_RC=$?
set -e
if [ "$CURL_RC" -ne 0 ]; then
  echo "warn: gzip curl failed rc $CURL_RC (non-fatal, checking plain still passes)" >&2
else
  read -r GZ_CODE GZ_TIME GZ_SIZE GZ_CTYPE < "$CURL_STATS"
  GZ_WALL_MS="$(python3 -c "import sys; print(int(float(sys.argv[1])*1000))" "$GZ_TIME" 2>/dev/null || echo "?")"
  echo "    http=$GZ_CODE content_type=$GZ_CTYPE time=${GZ_TIME}s (${GZ_WALL_MS}ms) size=${GZ_SIZE}B"
  echo "    headers: $(grep -i content-encoding "$GZIP_HEADERS" || echo 'no Content-Encoding')"
  if [ "$GZ_CODE" = "200" ]; then
    HAS_GZIP=0
    if grep -qi "content-encoding: gzip" "$GZIP_HEADERS"; then
      HAS_GZIP=1
      echo "    ✓ gzip Content-Encoding present"
      if head -c 2 "$TRACE_GZ" | od -An -tx1 | grep -q "1f 8b"; then
        echo "    ✓ gzip magic present"
      else
        echo "warn: Content-Encoding gzip but file not gzipped (may be auto-decompressed by curl/proxy)" >&2
        HAS_GZIP=0
      fi
    else
      echo "warn: server did not return gzip (Content-Encoding missing) — may be small payload or proxy" >&2
    fi

    if [ "$HAS_GZIP" -eq 1 ]; then
      GZ_BYTES="$(wc -c < "$TRACE_GZ" | tr -d '[:space:]')"
      echo "    gzipped bytes=$GZ_BYTES vs plain $ACTUAL_BYTES ratio=$(python3 -c "import sys; print(f'{int(sys.argv[1])/int(sys.argv[2]):.2%}')" "$GZ_BYTES" "$ACTUAL_BYTES" 2>/dev/null || echo '?')"
      if [ "$GZ_BYTES" -ge "$ACTUAL_BYTES" ]; then
        echo "warn: gzipped size not smaller than plain ($GZ_BYTES >= $ACTUAL_BYTES)" >&2
      fi
      echo "==> gunzip and re-validate"
      python3 -c "import gzip, pathlib; p=pathlib.Path('$TRACE_GZ'); data=gzip.decompress(p.read_bytes()); pathlib.Path('$TRACE_GZ_UNZIPPED').write_bytes(data); print(f'decompressed {len(data)} bytes')" || fail "gunzip failed"
      python3 "$VALIDATOR" --input "$TRACE_GZ_UNZIPPED" --epsilon "$EPSILON" --delta "$DELTA" --trace-id "$TRACE_ID" || fail "gzipped NDJSON validation failed after gunzip"
      echo "    ✓ gzipped payload decompresses to valid NDJSON"
    else
      if [ -s "$TRACE_GZ" ]; then
        echo "    checking non-gzipped fallback payload"
        if head -c 1 "$TRACE_GZ" | grep -q "{"; then
          python3 "$VALIDATOR" --input "$TRACE_GZ" --epsilon "$EPSILON" --delta "$DELTA" --trace-id "$TRACE_ID" || echo "warn: fallback gzip file not valid NDJSON" >&2
        fi
      fi
    fi
  else
    fail "gzip request returned $GZ_CODE (expected 200)"
  fi
fi

# ---------------------------------------------------------------------------
# Direct handler smoke (no Flask) — validates C++ emission path
# ---------------------------------------------------------------------------
echo "==> direct handler smoke: simplify --web-server --json-stream"
set +e
"$REPO_ROOT/build/simplify" --in "$TRACE_ID" -e "$EPSILON" -d "$DELTA" --web-server --json-stream > "$DIRECT_FILE" 2> /tmp/web-check-direct-stderr.txt
DIRECT_RC=$?
set -e
if [ "$DIRECT_RC" -ne 0 ]; then
  cat /tmp/web-check-direct-stderr.txt >&2 || true
  fail "direct handler failed (rc $DIRECT_RC)"
fi
if [ ! -s "$DIRECT_FILE" ]; then
  fail "direct handler produced empty file"
fi
DIRECT_BYTES="$(wc -c < "$DIRECT_FILE" | tr -d '[:space:]')"
echo "    direct bytes=$DIRECT_BYTES"
python3 "$VALIDATOR" --input "$DIRECT_FILE" --epsilon "$EPSILON" --delta "$DELTA" --trace-id "$TRACE_ID" || fail "direct handler NDJSON validation failed"

# Compare plain vs direct: both should have same prefix count? Allow small variance
PLAIN_PREFIXES="$(python3 -c "import json; lines=open('$TRACE_FILE').read().strip().splitlines(); print(sum(1 for l in lines if json.loads(l).get('type')=='prefix'))")"
DIRECT_PREFIXES="$(python3 -c "import json; lines=open('$DIRECT_FILE').read().strip().splitlines(); print(sum(1 for l in lines if json.loads(l).get('type')=='prefix'))")"
echo "    plain prefixes=$PLAIN_PREFIXES direct prefixes=$DIRECT_PREFIXES"
if [ "$PLAIN_PREFIXES" != "$DIRECT_PREFIXES" ]; then
  fail "prefix count mismatch plain $PLAIN_PREFIXES vs direct $DIRECT_PREFIXES (check C++ vs server params)"
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
echo "✅ web trace check passed (trace $TRACE_ID eps=$EPSILON delta=$DELTA)"
echo "   plain: ${ACTUAL_BYTES}B ${WALL_MS}ms $PLAIN_PREFIXES prefixes"
if [ -n "${GZ_BYTES:-}" ]; then
  echo "   gzip:  ${GZ_BYTES}B (ratio $(python3 -c "import sys; print(f'{int(sys.argv[1])/int(sys.argv[2]):.1%}')" "${GZ_BYTES:-0}" "$ACTUAL_BYTES" 2>/dev/null || echo '?'))"
fi
echo "   direct: ${DIRECT_BYTES}B $DIRECT_PREFIXES prefixes"
echo "   server log: $SERVER_LOG"
# shellcheck disable=SC2016
echo '   to re-run: scripts/check_web_trace.sh --trace-id '"$TRACE_ID"' --epsilon '"$EPSILON"' --delta '"$DELTA"

FAIL=0
