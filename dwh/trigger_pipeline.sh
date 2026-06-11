#!/usr/bin/env bash
#
# Trigger the DWH pipeline on Google Cloud Build.
#
# Usage:
#   ./trigger_pipeline.sh [--suppress-datasets id1,id2,...]
#   ./trigger_pipeline.sh --postgres-only --bq-source-suffix _gene_id_migration --pg-table-suffix _gene_id_migration
#
# Prerequisites:
#   - gcloud CLI installed and authenticated
#   - Environment variables set (via dev_secrets or equivalent):
#       GCLOUD_PROJECT, GCLOUD_REGION, BQ_DATASET, BQ_LOCATION, GCLOUD_TMP_BUCKET,
#       PG_CONN_INTERNAL, ES_URL, ES_USERNAME, ES_PASSWORD
#

set -euo pipefail

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
SUPPRESS_DATASETS=""
RUN_DBT="true"
RUN_POSTGRES="true"
RUN_ES="true"
BQ_SOURCE_SUFFIX="${BQ_SOURCE_SUFFIX:-}"
PG_TABLE_SUFFIX="${PG_TABLE_SUFFIX:-}"
SKIP_MATERIALIZED_VIEWS="${SKIP_MATERIALIZED_VIEWS:-false}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --suppress-datasets)
            SUPPRESS_DATASETS="$2"
            shift 2
            ;;
        --postgres-only)
            RUN_DBT="false"
            RUN_POSTGRES="true"
            RUN_ES="false"
            shift
            ;;
        --bq-source-suffix)
            BQ_SOURCE_SUFFIX="$2"
            shift 2
            ;;
        --pg-table-suffix)
            PG_TABLE_SUFFIX="$2"
            shift 2
            ;;
        --skip-materialized-views)
            SKIP_MATERIALIZED_VIEWS="true"
            shift
            ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: $0 [--suppress-datasets id1,id2,...] [--postgres-only] [--bq-source-suffix suffix] [--pg-table-suffix suffix] [--skip-materialized-views]"
            exit 1
            ;;
    esac
done

if [[ -n "$PG_TABLE_SUFFIX" ]]; then
    # Dev/suffixed Postgres syncs must not refresh production materialized views.
    SKIP_MATERIALIZED_VIEWS="true"
fi

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
)

if [[ "$RUN_ES" == "true" ]]; then
    REQUIRED_VARS+=(
    ES_URL
    ES_USERNAME
    ES_PASSWORD
)
fi

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
    echo "Make sure to source your secrets first, e.g.:  dev_secrets"
    exit 1
fi

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
echo "  BQ Source Suffix:   ${BQ_SOURCE_SUFFIX:-<none>}"
echo "  PG Table Suffix:    ${PG_TABLE_SUFFIX:-<none>}"
echo "  BQ Location:        $BQ_LOCATION"
echo "  GCS Bucket:         $GCLOUD_TMP_BUCKET"
echo "  Suppress Datasets:  ${SUPPRESS_DATASETS:-<none>}"
echo "  Run dbt:            $RUN_DBT"
echo "  Run Postgres:       $RUN_POSTGRES"
echo "  Run Elasticsearch:  $RUN_ES"
echo "  Skip MVs:           $SKIP_MATERIALIZED_VIEWS"
echo "============================================"
echo ""

BUILD_ID=$(gcloud builds submit "$SCRIPT_DIR" \
    --project="$GCLOUD_PROJECT" \
    --region="$GCLOUD_REGION" \
    --config="$SCRIPT_DIR/cloudbuild.yaml" \
    --gcs-source-staging-dir="gs://$GCLOUD_TMP_BUCKET/cloudbuild-source" \
    --substitutions="\
_GCLOUD_PROJECT=$GCLOUD_PROJECT,\
_GCLOUD_REGION=$GCLOUD_REGION,\
_BQ_DATASET=$BQ_DATASET,\
_BQ_LOCATION=$BQ_LOCATION,\
_GCLOUD_TMP_BUCKET=$GCLOUD_TMP_BUCKET,\
_PG_CONN_INTERNAL=$PG_CONN_INTERNAL,\
_ES_URL=${ES_URL:-},\
_ES_USERNAME=${ES_USERNAME:-},\
_ES_PASSWORD=${ES_PASSWORD:-},\
_SUPPRESS_DATASETS=$SUPPRESS_DATASETS,\
_RUN_DBT=$RUN_DBT,\
_RUN_POSTGRES=$RUN_POSTGRES,\
_RUN_ES=$RUN_ES,\
_BQ_SOURCE_SUFFIX=$BQ_SOURCE_SUFFIX,\
_PG_TABLE_SUFFIX=$PG_TABLE_SUFFIX,\
_SKIP_MATERIALIZED_VIEWS=$SKIP_MATERIALIZED_VIEWS" \
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
