# GitHub Actions Workflows

This directory contains CI/CD workflows for the trajectory simplification project.

## Workflows

### correctness.yml - Correctness regression
**Triggers:** Push to main, pull requests to main, manual dispatch

Compares Fréchet distance and point counts for derived benchmark trajectories against the baseline commit: the PR base (`github.event.pull_request.base.sha`) on pull requests, or the previous commit (`github.event.before`) on pushes.

Runs 20 parallel matrix jobs (same `(ε, δ, size)` triples as the benchmark gate):

| label | ε | δ formula | δ | size | IDs | why |
| --- | --- | --- | --- | --- | --- | --- |
| `extra-coarse-e-d300-small` | 299 | 300/(1+ε) | 1 | small ~100 pts | 11..20 | Largest ε, tight corridor (corridor area constant) |
| `extra-coarse-e-d300-large` | 299 | 300/(1+ε) | 1 | large ~1000 pts | 21..30 | Same ε/δ, larger trajectory stresses length scaling |
| `extra-coarse-e-d1000-small` | 299 | 1000/(1+ε) | 3.333… | small ~100 pts | 11..20 | Same ε, looser corridor (1000 scale) |
| `extra-coarse-e-d1000-large` | 299 | 1000/(1+ε) | 3.333… | large ~1000 pts | 21..30 | Same ε, larger + looser |
| `coarse-e-d300-small` | 30 | 300/(1+ε) | 9.677… | small | 11..20 | Coarse ε |
| `coarse-e-d300-large` | 30 | 300/(1+ε) | 9.677… | large | 21..30 | Coarse ε, large |
| `coarse-e-d1000-small` | 30 | 1000/(1+ε) | 32.258… | small | 11..20 | Coarse ε, looser δ |
| `coarse-e-d1000-large` | 30 | 1000/(1+ε) | 32.258… | large | 21..30 | Coarse ε, large + looser |
| `mid-e-d300-small` | 5 | 300/(1+ε) | 50 | small | 11..20 | Mid ε, δ=50 (historical mid) |
| `mid-e-d300-large` | 5 | 300/(1+ε) | 50 | large | 21..30 | Mid ε, large |
| `mid-e-d1000-small` | 5 | 1000/(1+ε) | 166.666… | small | 11..20 | Mid ε, looser δ |
| `mid-e-d1000-large` | 5 | 1000/(1+ε) | 166.666… | large | 21..30 | Mid ε, large + looser |
| `fine-e-d300-small` | 0.5 | 300/(1+ε) | 200 | small | 11..20 | Fine ε (headless default) |
| `fine-e-d300-large` | 0.5 | 300/(1+ε) | 200 | large | 21..30 | Fine ε, large |
| `fine-e-d1000-small` | 0.5 | 1000/(1+ε) | 666.666… | small | 11..20 | Fine ε, looser δ |
| `fine-e-d1000-large` | 0.5 | 1000/(1+ε) | 666.666… | large | 21..30 | Fine ε, large + looser |
| `extra-fine-e-d300-small` | 0.1 | 300/(1+ε) | 272.727… | small | 11..20 | Extra-fine ε |
| `extra-fine-e-d300-large` | 0.1 | 300/(1+ε) | 272.727… | large | 21..30 | Extra-fine ε, large |
| `extra-fine-e-d1000-small` | 0.1 | 1000/(1+ε) | 909.090… | small | 11..20 | Extra-fine ε, looser δ |
| `extra-fine-e-d1000-large` | 0.1 | 1000/(1+ε) | 909.090… | large | 21..30 | Extra-fine ε, large + looser |

Epsilon tiers are `extra-coarse-e` 299, `coarse-e` 30, `mid-e` 5, `fine-e` 0.5, `extra-fine-e` 0.1. Delta is computed exactly as `300/(1+e)` or `1000/(1+e)` in the job (`python3 -c "numer/(1+eps)"` with `:.15g`), not hardcoded, so floating error cannot drift. The `d300` / `d1000` suffix and `small`/`large` suffix in the label make each of the 20 configurations self-describing. Size maps to deterministic datasets derived by `scripts/derive_benchmark_data.py`: `small` IDs 11..20 are first 100 points of original 1..10 (data/8 stays 5, avg 90.5), `large` IDs 21..30 are first 1000 points of original 1..10 (3×1000 + 7 smaller, avg ~695). Original data/1..10 sizes 5..5618 are not uniformly ~100 or ~1000, so the prefix derivation is required; see the script header for rationale and alternatives considered.

Tolerances (unchanged without evidence): `DIST_TOL=0.01`, `POINTS_TOL=0`.

Artifacts: `correctness-<label>` with `correctness.tsv` and `correctness.json`. Job summary tables show per-ID wins/fails. Delta is printed as `DELTA_NUMER/(1+EPSILON) = DELTA` in logs.

### benchmark.yml - Performance regression
**Triggers:** Push to main, pull requests to main, manual dispatch

Same 20 parallel `(ε, δ, size)` matrix jobs as correctness (table above). For each setting, runs `BENCH_RUNS=10` Release invocations of `SIMPLIFY_CORE_MS` per ID for each binary (small: 11..20 or large: 21..30, new vs the same baseline commit as correctness), computes per-ID mean and sample stddev via Welford's online algorithm, and enforces high-confidence gates (one-sided 95% CI with Welch t, fail only when confident new is worse):

1. **Mean gate:** over IDs with `orig_ms ≥ MIN_BENCH_MS` (1 ms), `mean(new) ≤ mean(orig) × MEAN_LIMIT` with `MEAN_LIMIT=1.05` (cannot worsen by more than 1.05×) but gated as `(mu_new - 1.05·mu_orig) - t_{0.95}·SE_mean > 0` where `SE_mean` propagates per-ID `σ/√n` via Welch-Satterthwaite. Replaces a zero-tolerance `mean(new) < mean(orig)` check that failed `mid` on ~0.2% noise with `gated_n=1` and the prior naive 1.05× on 5 runs that still failed spuriously 80/82.
2. **Per-ID gate:** for those same gated IDs, `new_ms ≤ orig_ms * 1.20` (was 1.50) but now high-confidence: `SE = sqrt(σ_new²/n + (1.20·σ_orig)²/n)`, Welch `df`, `t_{0.95}`, fail `FAIL_SLOW` iff `(mu_new - 1.20·mu_orig) - t·SE > 0`.

IDs with `orig_ms` strictly below 1 ms are reported as `SKIP_FLOOR` and excluded from both gates; everything at or above 1 ms is gated. Before the confidence gate the 5-run naive means failed 80/82 matrix runs spuriously on 5-run noise; the 10-run Welch gating eliminates those false positives while still catching true ~20%+ regressions.

Local `SIMPLIFY_CORE_MS` means can swing under host load. When comparing commits locally, interleave baseline and candidate runs and inspect per-ID minima alongside the ten-run Welford means/stddevs used by CI; the benchmark JSON/TSV now include `orig_std`/`new_std` per ID.

Gated averages intentionally omit `--time` so timing stays comparable to older binaries. After averages, the new binary runs once per ID with `--time` and must exit 0 with at least one `TIMER_MS` line (blank ops is a failure). Phase counters (`hull_Gi`, `find_F`, `intersect`, `boundary_P`, …) land in the TSV `ops` column and nested `ops` objects in JSON.

Artifacts: `benchmark-<label>` with `benchmark.tsv` / `benchmark.json`. TSV header is `id e d orig_ms orig_std new_ms new_std ratio gated status ops`; JSON `cases` include `orig_std`/`new_std` (sample stddev, Welford) alongside `orig_ms`/`new_ms`. Job summaries show `orig_ms ± std` / `new_ms ± std` and mark `FAIL_SLOW` only when the Welch 95% lower bound exceeds the limit. Each job logs `Computed DELTA=… from NUMER/(1+EPSILON)` and `Synced IDs: …`.

Runtime: 20 jobs each build both binaries and run 10 IDs × (10 orig + 10 new + 1 new --time) = 210 simplify invocations per job. On a local Release build the full 20×10 correctness sweep (~200 Fréchet runs) and benchmark sweep (4,000 sampled invocations plus 200 ops runs, 4,200 total) each complete in under a few minutes; measured total for the entire matrix (both workflows sequentially, `BENCH_RUNS=10`) is reported in the PR that introduced the matrix. The prior `BENCH_RUNS=5` naive-mean gate failed 80/82 runs spuriously; the 10-run high-confidence gate removes that noise.

### gui-build.yml - Qt GUI build
**Triggers:** Push to main, pull requests to main, manual dispatch

Builds the `simplify_with_gui` target with the default `BUILD_GUI=ON`
using distro Qt packages (`qt6-base-dev`, `qt6-svg-dev`) plus pinned
CGAL 6.x headers fetched from the CGAL release tarball (cached): Ubuntu
24.04 ships CGAL 5.6, which is Qt5-only, while the GUI target needs
`CGAL::CGAL_Qt6`. The correctness and benchmark workflows build with
`-DBUILD_GUI=OFF`, so without this job a Qt-only breakage (e.g. an
identifier colliding with Qt's `emit` macro) compiles clean in CI.

### deploy.yml - Cloud Run deploy
**Triggers:** Push to main, manual dispatch

Builds the repo `Dockerfile` with Cloud Build (`gcloud run deploy --source .`) and publishes the web viewer to Cloud Run.

Cloud Build can flake on network-bound Julia install / `Pkg.add` (opaque "Building Container ... failed" after README-only #21). The `Dockerfile` retries `curl`/juliaup install and `Pkg.add` with loud failure logs; image behavior is unchanged.

Authentication uses Workload Identity Federation (OIDC). No repository secret is required. The workflow requests `id-token: write` and impersonates the deploy service account through a GitHub-restricted identity pool provider.

Defaults (override with repository variables):

| Variable | Default |
| --- | --- |
| `GCP_PROJECT_ID` | `project-ec366840-6857-446e-852` |
| `GCP_REGION` | `asia-east2` |
| `CLOUD_RUN_SERVICE` | `simplify-viewer` |
| `GCP_WORKLOAD_IDENTITY_PROVIDER` | `projects/522405269791/locations/global/workloadIdentityPools/github-actions/providers/github` |
| `GCP_SERVICE_ACCOUNT` | `github-actions-deploy@project-ec366840-6857-446e-852.iam.gserviceaccount.com` |

GCP setup already applied for this repo (`yeungsinchun/Simplification-of-Trajectory-Streams`):

- Workload Identity Pool `github-actions` with OIDC provider `github` (attribute condition locks to this repository)
- Service account `github-actions-deploy@…` with Cloud Run Admin, Cloud Build Editor, Service Account User, Storage Admin, and Artifact Registry Admin
- Pool principal bound as `roles/iam.workloadIdentityUser` on that service account

Service flags match the former local `deploy.sh`: 4 GiB RAM, 2 CPU, 300s timeout, max 10 instances, `--no-cpu-throttling` (needed so background Julia Fréchet work keeps CPU after `/api/frechet` returns), `--allow-unauthenticated`.

## Local deploy

Keep a machine-local `deploy.sh` (gitignored) or run the same `gcloud run deploy` command from `deploy.yml`.
