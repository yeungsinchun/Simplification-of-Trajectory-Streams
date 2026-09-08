# GitHub Actions Workflows

This directory contains CI/CD workflows for the trajectory simplification project.

## Workflows

### correctness.yml - Correctness regression
**Triggers:** Push to main, pull requests to main, manual dispatch

Compares Fréchet distance and point counts for datasets 1..10 against the previous commit.

Runs four parallel matrix jobs (same `(ε, δ)` pairs as the benchmark gate):

| label | ε | δ | why |
| --- | --- | --- | --- |
| `coarse-e` | 299 | 1 | Historical CI pair (large match tolerance) |
| `fine-e` | 0.5 | 300 | Tight ε (headless CLI default) with large δ (corridor constant from `scripts/benchmark_e.py`) |
| `mid` | 5 | 50 | Mid-range; δ ≈ 300/(1+ε) corridor scaling from `scripts/benchmark_e.py` |
| `large-d` | 1 | 1000 | Small ε with a very large δ (stresses the time/corridor axis differently from `fine-e`) |

Tolerances (unchanged without evidence): `DIST_TOL=0.01`, `POINTS_TOL=0`.

Artifacts: `correctness-<label>-e…-d…` with `correctness.tsv` and `correctness.json`. Job summary tables show per-ID wins/fails.

### benchmark.yml - Performance regression
**Triggers:** Push to main, pull requests to main, manual dispatch

Same four parallel `(ε, δ)` matrix jobs as correctness. For each setting, averages `BENCH_RUNS=5` Release runs of `SIMPLIFY_CORE_MS` on IDs 1..10 (new vs previous commit) and enforces:

1. **Mean gate:** over IDs with `orig_ms ≥ MIN_BENCH_MS` (20 ms), `mean(new) ≤ mean(orig) × MEAN_LIMIT` with `MEAN_LIMIT=1.05` (cannot worsen by more than 1.05×). Replaces a zero-tolerance `mean(new) < mean(orig)` check that failed `mid` on ~0.2% noise with `gated_n=1`.
2. **Per-ID gate:** for those same gated IDs, `new_ms ≤ orig_ms * 1.20` (was 1.50).

IDs below the 20 ms floor are reported as `SKIP_FLOOR` and excluded from both gates so wall-clock noise does not fail the job.

Gated averages intentionally omit `--time` so timing stays comparable to older binaries. After averages, the new binary runs once per ID with `--time` and must exit 0 with at least one `TIMER_MS` line (blank ops is a failure). Phase counters (`hull_Gi`, `find_F`, `intersect`, `boundary_P`, …) land in the TSV `ops` column and nested `ops` objects in JSON.

Artifacts: `benchmark-<label>-e…-d…` with `benchmark.tsv` / `benchmark.json`. Job summaries highlight wins and regressions per setting.

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
