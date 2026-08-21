#!/bin/bash
set -e

# Cloud Run deployment script for trajectory simplification visualizer

PROJECT_ID="${GOOGLE_CLOUD_PROJECT:-project-ec366840-6857-446e-852}"
REGION="${REGION:-asia-east2}"
SERVICE_NAME="simplify-viewer"

echo "==> Deploying to Google Cloud Run"
echo "    Project: $PROJECT_ID"
echo "    Region: $REGION"
echo "    Service: $SERVICE_NAME"
echo ""

# Build and deploy using Cloud Build
# --no-cpu-throttling: CPU is always allocated (not gated to request windows).
#   Required because /api/frechet returns immediately and does the heavy
#   Julia work in a background daemon thread. Without this, the Julia
#   process never gets CPU time after the request completes.
# --memory 4Gi: trace 1 (588 pts) + CGAL + Julia FrechetDist can exceed 2 GiB.
gcloud run deploy "$SERVICE_NAME" \
  --source . \
  --platform managed \
  --region "$REGION" \
  --project "$PROJECT_ID" \
  --allow-unauthenticated \
  --memory 4Gi \
  --cpu 2 \
  --timeout 300 \
  --max-instances 10 \
  --no-cpu-throttling

echo ""
echo "==> Deployment complete!"
echo "    URL: https://$SERVICE_NAME-$(gcloud run services describe $SERVICE_NAME --region $REGION --format='value(status.url)' | cut -d/ -f3)"
