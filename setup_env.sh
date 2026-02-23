#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────
#  Setup script for Breast Cancer scRNA-seq TME Analysis
# ──────────────────────────────────────────────────────────
set -euo pipefail

VENV_DIR=".venv"
KERNEL_NAME="breast-cancer-scrna"
PYTHON="${PYTHON:-python3}"

# ── Helpers ───────────────────────────────────────────────
info()  { printf "\033[1;34m[INFO]\033[0m  %s\n" "$*"; }
ok()    { printf "\033[1;32m[OK]\033[0m    %s\n" "$*"; }
warn()  { printf "\033[1;33m[WARN]\033[0m  %s\n" "$*"; }
error() { printf "\033[1;31m[ERROR]\033[0m %s\n" "$*"; exit 1; }

# ── Pre-flight checks ────────────────────────────────────
command -v "$PYTHON" >/dev/null 2>&1 || error "Python not found. Install Python 3.10+ first."

PY_VERSION=$("$PYTHON" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
PY_MAJOR=$("$PYTHON" -c 'import sys; print(sys.version_info.major)')
PY_MINOR=$("$PYTHON" -c 'import sys; print(sys.version_info.minor)')

if [[ "$PY_MAJOR" -lt 3 ]] || [[ "$PY_MAJOR" -eq 3 && "$PY_MINOR" -lt 10 ]]; then
    error "Python 3.10+ required (found ${PY_VERSION})"
fi

info "Using Python ${PY_VERSION} ($(command -v "$PYTHON"))"

# ── Create virtual environment ────────────────────────────
if [[ -d "$VENV_DIR" ]]; then
    warn "Virtual environment already exists at ${VENV_DIR}/"
    read -rp "  Recreate it? [y/N] " answer
    if [[ "${answer,,}" == "y" ]]; then
        info "Removing existing virtual environment..."
        rm -rf "$VENV_DIR"
    else
        info "Reusing existing virtual environment."
    fi
fi

if [[ ! -d "$VENV_DIR" ]]; then
    info "Creating virtual environment in ${VENV_DIR}/..."
    "$PYTHON" -m venv "$VENV_DIR"
    ok "Virtual environment created."
fi

# ── Activate ──────────────────────────────────────────────
# shellcheck disable=SC1091
source "${VENV_DIR}/bin/activate"
info "Activated virtual environment."

# ── Upgrade pip ───────────────────────────────────────────
info "Upgrading pip..."
pip install --upgrade pip --quiet

# ── Install dependencies ─────────────────────────────────
info "Installing dependencies from requirements.txt..."
pip install -r requirements.txt --quiet

# scikit-misc is needed for scanpy's seurat_v3 HVG method
info "Installing scikit-misc (required by seurat_v3 HVG selection)..."
pip install scikit-misc --quiet

ok "All packages installed."

# ── Register Jupyter kernel ───────────────────────────────
info "Registering Jupyter kernel '${KERNEL_NAME}'..."
python -m ipykernel install --user \
    --name "$KERNEL_NAME" \
    --display-name "Breast Cancer scRNA-seq (Python ${PY_VERSION})" \
    --quiet 2>/dev/null || \
python -m ipykernel install --user \
    --name "$KERNEL_NAME" \
    --display-name "Breast Cancer scRNA-seq (Python ${PY_VERSION})"

ok "Jupyter kernel registered."

# ── Create data directories ──────────────────────────────
mkdir -p data/raw data/processed results reports
ok "Directory structure verified."

# ── Summary ───────────────────────────────────────────────
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
ok "Environment setup complete!"
echo ""
echo "  To activate the environment:"
echo "    source ${VENV_DIR}/bin/activate"
echo ""
echo "  To run the notebook:"
echo "    jupyter notebook notebooks/01_breast_cancer_scrna_tme_analysis.ipynb"
echo ""
echo "  To run scripts directly:"
echo "    python scripts/01_download_and_qc.py"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
