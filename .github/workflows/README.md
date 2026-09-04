# GitHub Actions Workflows

This directory contains CI/CD workflows for the trajectory simplification project.

## Workflows

### correctness.yml - Correctness regression
**Triggers:** Push to main, pull requests to main, manual dispatch

Compares Fréchet distance and point counts for datasets 1..10 against the previous commit.

### benchmark.yml - Comprehensive Benchmarks
**Triggers:** Push to main, pull requests, daily schedule, manual dispatch

Builds in Release mode, sweeps epsilon values, and stores historical results.

### deploy.yml - Cloud Run deploy
**Triggers:** Push to main, manual dispatch

Builds the repo `Dockerfile` with Cloud Build (`gcloud run deploy --source .`) and publishes the web viewer to Cloud Run.

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

## Local Testing

To run the same Frechet, point-count, and performance checks locally, use
`scripts/local_ci_gate.py`. See that script's module docstring and `--help`.

## Local deploy

Keep a machine-local `deploy.sh` (gitignored) or run the same `gcloud run deploy` command from `deploy.yml`.
