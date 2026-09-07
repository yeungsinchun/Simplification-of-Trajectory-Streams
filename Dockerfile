FROM ubuntu:22.04 AS builder

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libcgal-dev \
    qt6-base-dev \
    git \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# Headless simplify sources
COPY simplify.cpp simplify_geometry.h simplify_io.h timer.h CMakeLists.txt ./

# Pin traj-compression to the gitlink SHA. The Docker context excludes .git
# (and Cloud Run uploads exclude the submodule), so baseline sources are
# fetched here. BUILD_GUI=OFF skips the Qt viewer; dots still builds via
# Qt6 Core from qt6-base-dev.
RUN git clone https://github.com/yeungsinchun/traj-compression.git traj-compression && \
    git -C traj-compression checkout ce40df79a77db8eb4bb31d7316a187c07b2ad921

RUN cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_GUI=OFF && \
    cmake --build build --target simplify dots dp squish -j"$(nproc)"

# Final stage: Python runtime on Ubuntu
FROM ubuntu:22.04

# Install Python and runtime dependencies.
# CGAL is header-only at build time; simplify needs GMP/MPFR at runtime.
# dots needs Qt6 Core (libQt6Core).
RUN apt-get update && apt-get install -y \
    python3.11 \
    python3-pip \
    libgmp10 \
    libmpfr6 \
    libqt6core6 \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install Julia via the official install script. It auto-detects the current
# platform architecture (x86_64 vs aarch64), downloads the right tarball via
# juliaup, and exposes the `julia` binary on PATH. The default install root
# is /root/.juliaup when running as root, which we add to PATH for later
# steps. The Ubuntu 22.04 apt repo does not ship a usable `julia` package,
# so this is the canonical install path.
#
# Cloud Build has flaked here with opaque "Building Container ... failed"
# (often curl/network during juliaup or later Pkg.add). Retry downloads and
# fail loudly so the next Cloud Build log shows the real error.
ENV JULIAUP_INSTALL=/root/.juliaup
ENV PATH="$JULIAUP_INSTALL/bin:$PATH"
RUN set -eux; \
    ok=0; \
    for attempt in 1 2 3 4 5; do \
      echo "juliaup install: download script attempt ${attempt}/5"; \
      if curl --retry 5 --retry-delay 5 --retry-connrefused -fsSL \
           https://install.julialang.org -o /tmp/julia-install.sh; then \
        ok=1; \
        break; \
      fi; \
      echo "juliaup install: curl failed on attempt ${attempt}; sleeping..."; \
      sleep $((attempt * 5)); \
    done; \
    if [ "$ok" -ne 1 ]; then \
      echo "ERROR: failed to download https://install.julialang.org after 5 attempts" >&2; \
      exit 1; \
    fi; \
    ok=0; \
    for attempt in 1 2 3 4 5; do \
      echo "juliaup install: run installer attempt ${attempt}/5"; \
      if sh /tmp/julia-install.sh --yes; then \
        ok=1; \
        break; \
      fi; \
      echo "juliaup install: installer failed on attempt ${attempt}; sleeping..."; \
      sleep $((attempt * 5)); \
    done; \
    rm -f /tmp/julia-install.sh; \
    if [ "$ok" -ne 1 ]; then \
      echo "ERROR: juliaup installer failed after 5 attempts" >&2; \
      exit 1; \
    fi; \
    command -v julia; \
    julia --version

# Sanity-check that the julia binary actually runs in the final image
RUN julia --version

# Make python3.11 the default
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.11 1 && \
    update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1

WORKDIR /app

# Copy compiled binaries from builder (web compare pane needs all four)
COPY --from=builder /build/build/simplify /app/build/simplify
COPY --from=builder /build/build/dots /app/build/dots
COPY --from=builder /build/build/dp /app/build/dp
COPY --from=builder /build/build/squish /app/build/squish

# Copy web application files
COPY web/ /app/web/
COPY data/ /app/data/
COPY scripts/ /app/scripts/

# Install Python dependencies
RUN pip install --no-cache-dir -r /app/web/requirements.txt

# Pre-install the FrechetDist and ArgParse Julia packages so the first
# /api/frechet request doesn't pay the full Pkg.add cost (and so the image is
# self-contained). ArgParse is required by scripts/frechet.jl for --batch and
# --raw flag parsing. Retry Pkg.add - registry fetch flakes in Cloud Build.
RUN set -eux; \
    ok=0; \
    for attempt in 1 2 3 4 5; do \
      echo "julia Pkg.add: attempt ${attempt}/5"; \
      if julia -e 'using Pkg; Pkg.add(["FrechetDist", "ArgParse"])'; then \
        ok=1; \
        break; \
      fi; \
      echo "julia Pkg.add: failed on attempt ${attempt}; sleeping..."; \
      sleep $((attempt * 5)); \
    done; \
    if [ "$ok" -ne 1 ]; then \
      echo "ERROR: julia Pkg.add([FrechetDist, ArgParse]) failed after 5 attempts" >&2; \
      exit 1; \
    fi; \
    julia -e 'using FrechetDist, ArgParse; println("julia packages ok")'

# Set environment variables.
# Image default PORT=5050 avoids macOS AirPlay on 5000 when published locally.
# Local `python web/server.py` defaults to 5051 when PORT is unset; Cloud Run
# still injects PORT (typically 8080) at runtime.
ENV PORT=5050
ENV PYTHONUNBUFFERED=1

EXPOSE 5050

# Run the Flask server
CMD ["python", "/app/web/server.py"]
