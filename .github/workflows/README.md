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

Defaults (override with repository variables):

| Variable | Default |
| --- | --- |
| `GCP_PROJECT_ID` | `project-ec366840-6857-446e-852` |
| `GCP_REGION` | `asia-east2` |
| `CLOUD_RUN_SERVICE` | `simplify-viewer` |

Required repository secret:

- `GCP_SA_KEY` - JSON key for a service account that can deploy Cloud Run from source (Cloud Run Admin, Cloud Build Editor, Service Account User, and Storage access for source upload). Enable the Cloud Run, Cloud Build, and Artifact Registry APIs in the project.

Service flags match the former local `deploy.sh`: 4 GiB RAM, 2 CPU, 300s timeout, max 10 instances, `--no-cpu-throttling` (needed so background Julia Fréchet work keeps CPU after `/api/frechet` returns), `--allow-unauthenticated`.

## Local deploy

Keep a machine-local `deploy.sh` (gitignored) or run the same `gcloud run deploy` command from `deploy.yml`.
