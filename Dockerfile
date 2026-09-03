FROM ubuntu:22.04 AS builder

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libcgal-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy source code
WORKDIR /build
COPY simplify.cpp simplify_geometry.h simplify_io.h timer.h CMakeLists.txt ./

# Build the headless simplify binary only. BUILD_GUI=OFF skips the Qt viewer
# and the DOTS baseline so the resulting image has no Qt6 dependency.
RUN cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_GUI=OFF && \
    cmake --build build --target simplify

# Final stage: Python runtime on Ubuntu
FROM ubuntu:22.04

# Install Python and runtime dependencies
# Note: CGAL (header-only) and Julia ship their own dependencies.
# The simplify binary only needs GMP/MPFR; Boost is not needed at runtime.
RUN apt-get update && apt-get install -y \
    python3.11 \
    python3-pip \
    libgmp10 \
    libmpfr6 \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install Julia via the official install script. It auto-detects the current
# platform architecture (x86_64 vs aarch64), downloads the right tarball via
# juliaup, and exposes the `julia` binary on PATH. The default install root
# is /root/.juliaup when running as root, which we add to PATH for later
# steps. The Ubuntu 22.04 apt repo does not ship a usable `julia` package,
# so this is the canonical install path.
ENV JULIAUP_INSTALL=/root/.juliaup
ENV PATH="$JULIAUP_INSTALL/bin:$PATH"
RUN curl -fsSL https://install.julialang.org -o /tmp/julia-install.sh \
    && sh /tmp/julia-install.sh --yes \
    && rm /tmp/julia-install.sh

# Sanity-check that the julia binary actually runs in the final image
RUN julia --version

# Make python3.11 the default
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.11 1 && \
    update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1

WORKDIR /app

# Copy the compiled binary from builder
COPY --from=builder /build/build/simplify /app/build/simplify

# Copy web application files
COPY web/ /app/web/
COPY data/ /app/data/
COPY scripts/ /app/scripts/

# Install Python dependencies
RUN pip install --no-cache-dir -r /app/web/requirements.txt

# Pre-install the FrechetDist and ArgParse Julia packages so the first
# /api/frechet request doesn't pay the full Pkg.add cost (and so the image is
# self-contained). ArgParse is required by scripts/frechet.jl for --batch and
# --raw flag parsing.
RUN julia -e 'using Pkg; Pkg.add(["FrechetDist", "ArgParse"])'

# Set environment variables
# Default to 5050 to match local Flask (`python server.py`) and avoid macOS
# AirPlay on 5000. Cloud Run still injects PORT (typically 8080) at runtime.
ENV PORT=5050
ENV PYTHONUNBUFFERED=1

EXPOSE 5050

# Run the Flask server
CMD ["python", "/app/web/server.py"]
