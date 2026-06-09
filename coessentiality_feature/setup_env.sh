#!/usr/bin/env bash
# Sets up a Python virtual environment and installs all dependencies.
# Run once from the coessentiality_feature/ directory:
#   bash setup_env.sh
#
# Activate later with:
#   source .venv/bin/activate

set -euo pipefail

PYTHON=${PYTHON:-python3}
VENV_DIR=".venv"
REQUIREMENTS="requirements.txt"

# ── Python version check (need 3.10+) ─────────────────────────────────────
PY_VERSION=$("$PYTHON" -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
PY_MAJOR=$(echo "$PY_VERSION" | cut -d. -f1)
PY_MINOR=$(echo "$PY_VERSION" | cut -d. -f2)

if [ "$PY_MAJOR" -lt 3 ] || { [ "$PY_MAJOR" -eq 3 ] && [ "$PY_MINOR" -lt 10 ]; }; then
    echo "ERROR: Python 3.10+ is required. Found $PY_VERSION." >&2
    exit 1
fi
echo "Python $PY_VERSION detected."

# ── Create virtual environment ─────────────────────────────────────────────
if [ -d "$VENV_DIR" ]; then
    echo "Virtual environment already exists at $VENV_DIR — skipping creation."
else
    echo "Creating virtual environment at $VENV_DIR ..."
    "$PYTHON" -m venv "$VENV_DIR"
fi

# ── Activate and install ───────────────────────────────────────────────────
# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"

echo "Upgrading pip ..."
pip install --quiet --upgrade pip

echo "Installing dependencies from $REQUIREMENTS ..."
pip install --quiet -r "$REQUIREMENTS"

echo ""
echo "Done. Activate the environment with:"
echo "  source $VENV_DIR/bin/activate"
echo ""
echo "Then run the app with:"
echo "  python app/coessentiality_feature_pc.py"
