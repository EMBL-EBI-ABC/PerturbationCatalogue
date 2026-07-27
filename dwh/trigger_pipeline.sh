#!/usr/bin/env bash
#
# Trigger the DWH pipeline on Google Cloud Build.
#
# Usage:
#   ./trigger_pipeline.sh [--suppress-datasets id1,id2,...]
#
# Prerequisites:
#   - gcloud CLI installed and authenticated
#   - Environment variables set (via pc_secrets dev or equivalent):
#       GCLOUD_PROJECT, GCLOUD_REGION, BQ_DATASET, BQ_LOCATION, GCLOUD_TMP_BUCKET,
#       PG_CONN_INTERNAL, ES_URL, ES_USERNAME, ES_PASSWORD
#       Optional: ES_INDEX_SET (defaults to empty), OPENTARGETS_RELEASE (defaults to 26.03)
#

set -euo pipefail

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
SUPPRESS_DATASETS=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --suppress-datasets)
            SUPPRESS_DATASETS="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: $0 [--suppress-datasets id1,id2,...]"
            exit 1
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Validate environment
# ---------------------------------------------------------------------------
REQUIRED_VARS=(
    GCLOUD_PROJECT
    GCLOUD_REGION
    BQ_DATASET
    BQ_LOCATION
    GCLOUD_TMP_BUCKET
    PG_CONN_INTERNAL
    ES_URL
    ES_USERNAME
    ES_PASSWORD
)

missing=()
for var in "${REQUIRED_VARS[@]}"; do
    if [[ -z "${!var:-}" ]]; then
        missing+=("$var")
    fi
done

if [[ ${#missing[@]} -gt 0 ]]; then
    echo "ERROR: The following required environment variables are not set:"
    printf '  %s\n' "${missing[@]}"
    echo ""
    echo "Make sure to source your secrets first, e.g.:  pc_secrets dev"
    exit 1
fi

if [[ "$GCLOUD_PROJECT" == *"prod"* ]]; then
    echo "ERROR: Refusing to run the pipeline in a production environment: $GCLOUD_PROJECT" >&2
    exit 1
fi

OPENTARGETS_RELEASE="${OPENTARGETS_RELEASE:-26.03}"
ES_INDEX_SET="${ES_INDEX_SET:-}"

# ---------------------------------------------------------------------------
# Resolve paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# ---------------------------------------------------------------------------
# Copy files needed from outside dwh/ into the build context
# ---------------------------------------------------------------------------
mkdir -p "$SCRIPT_DIR/be"
cp "$REPO_ROOT/be/dataset_metadata.json" "$SCRIPT_DIR/be/dataset_metadata.json"
trap 'rm -rf "$SCRIPT_DIR/be"' EXIT

# ---------------------------------------------------------------------------
# Submit Cloud Build
# ---------------------------------------------------------------------------
echo "============================================"
echo " DWH Pipeline — Cloud Build"
echo "============================================"
echo "  Project:            $GCLOUD_PROJECT"
echo "  Region:             $GCLOUD_REGION"
echo "  BQ Dataset:         $BQ_DATASET"
echo "  BQ Reference:       $BQ_DATASET.opentargets_targets"
echo "  OT Release:         $OPENTARGETS_RELEASE"
echo "  BQ Location:        $BQ_LOCATION"
echo "  GCS Bucket:         $GCLOUD_TMP_BUCKET"
echo "  ES Index Set:       ${ES_INDEX_SET:-<default>}"
echo "  Suppress Datasets:  ${SUPPRESS_DATASETS:-<none>}"
echo "============================================"
echo ""

BUILD_ID=$(gcloud builds submit "$SCRIPT_DIR" \
    --project="$GCLOUD_PROJECT" \
    --region="$GCLOUD_REGION" \
    --config="$SCRIPT_DIR/cloudbuild.yaml" \
    --gcs-source-staging-dir="gs://$GCLOUD_TMP_BUCKET/cloudbuild-source" \
    --substitutions="^|^\
_GCLOUD_PROJECT=$GCLOUD_PROJECT|\
_GCLOUD_REGION=$GCLOUD_REGION|\
_BQ_DATASET=$BQ_DATASET|\
_OPENTARGETS_RELEASE=$OPENTARGETS_RELEASE|\
_BQ_LOCATION=$BQ_LOCATION|\
_GCLOUD_TMP_BUCKET=$GCLOUD_TMP_BUCKET|\
_PG_CONN_INTERNAL=$PG_CONN_INTERNAL|\
_ES_URL=$ES_URL|\
_ES_USERNAME=$ES_USERNAME|\
_ES_PASSWORD=$ES_PASSWORD|\
_ES_INDEX_SET=$ES_INDEX_SET|\
_SUPPRESS_DATASETS=$SUPPRESS_DATASETS" \
    --async \
    --format='value(id)')

echo ""
echo "Build submitted: $BUILD_ID"
echo ""
echo "Streaming logs (safe to interrupt — build continues in the cloud)..."
echo "To re-attach later:  gcloud builds log --stream $BUILD_ID --region=$GCLOUD_REGION --project=$GCLOUD_PROJECT"
echo ""

# Stream logs. If interrupted (e.g. laptop sleep), the build continues.
gcloud builds log --stream "$BUILD_ID" --region="$GCLOUD_REGION" --project="$GCLOUD_PROJECT" || true

# Check final status
STATUS=$(gcloud builds describe "$BUILD_ID" \
    --project="$GCLOUD_PROJECT" \
    --region="$GCLOUD_REGION" \
    --format='value(status)')

echo ""
echo "============================================"
echo " Build finished with status: $STATUS"
echo "============================================"

if [[ "$STATUS" != "SUCCESS" ]]; then
    echo "Build did not succeed. Check logs:"
    echo "  gcloud builds log $BUILD_ID --region=$GCLOUD_REGION --project=$GCLOUD_PROJECT"
    exit 1
fi
