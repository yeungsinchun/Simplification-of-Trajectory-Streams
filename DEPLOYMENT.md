# Cloud Run Deployment Guide

This document explains how to deploy the trajectory simplification visualizer to Google Cloud Run.

## Prerequisites

1. Install [Google Cloud SDK](https://cloud.google.com/sdk/docs/install)
2. Authenticate: `gcloud auth login`
3. Set your project: `gcloud config set project YOUR_PROJECT_ID`

## Quick Deploy

```bash
# Set your project ID
export GOOGLE_CLOUD_PROJECT="your-project-id-here"

# Optional: change region (defaults to us-central1)
export REGION="us-central1"

# Deploy
./deploy.sh
```

The script will:
- Build the Docker container using Cloud Build
- Deploy to Cloud Run with 2GB memory, 2 CPU
- Set up auto-scaling (max 10 instances)
- Allow unauthenticated access
- Return a public HTTPS URL

## Manual Deploy

If you prefer to deploy manually:

```bash
gcloud run deploy simplify-viewer \
  --source . \
  --platform managed \
  --region us-central1 \
  --allow-unauthenticated \
  --memory 2Gi \
  --cpu 2 \
  --timeout 300
```

## Configuration

The service uses these environment variables:
- `PORT`: Set automatically by Cloud Run (8080)
- `HOST`: Set to `0.0.0.0` for Cloud Run (reads from env)

## Cost Estimate

Cloud Run pricing (as of 2026):
- Free tier: 2M requests/month
- CPU: ~$0.00002400/vCPU-second
- Memory: ~$0.00000250/GiB-second
- Scales to zero when idle (no charges)

Typical cost for moderate use: **$5-20/month**

## Local Testing

Test the Docker container locally before deploying:

```bash
# Build the image
docker build -t simplify-viewer .

# Run locally
docker run -p 8080:8080 -e HOST=0.0.0.0 -e PORT=8080 simplify-viewer

# Open http://localhost:8080
```

## Troubleshooting

**Build fails**: Check that all source files are present and CGAL dependencies install correctly.

**Service timeout**: Increase with `--timeout 600` (max 3600 seconds).

**Out of memory**: Increase with `--memory 4Gi` or `--memory 8Gi`.

**View logs**: 
```bash
gcloud run services logs read simplify-viewer --region us-central1
```
