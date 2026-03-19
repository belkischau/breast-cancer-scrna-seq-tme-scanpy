#!/usr/bin/env python3
"""
01_download_and_qc.py
=====================
Download a public breast cancer scRNA-seq dataset and perform quality control.

Data source: GSE176078 — Wu et al., Nature Genetics 2021
"A single-cell and spatially resolved atlas of human breast cancers"

Two run modes:
  --mode full        Download / load all cells
  --mode quick-test  Subset to ~5,000 cells for fast demo (< 5 min)

If the GEO download fails (firewall, timeout, etc.), a high-fidelity synthetic
dataset is generated automatically so the pipeline can still be demonstrated.

Usage:
    python scripts/01_download_and_qc.py --mode quick-test
    python scripts/01_download_and_qc.py --mode full
"""

from __future__ import annotations

import argparse
import gzip
import io
import os
import shutil
import sys
import tarfile
import tempfile
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import nbinom

warnings.filterwarnings("ignore", category=FutureWarning)

# ──────────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"

QC_H5AD = PROCESSED_DIR / "breast_cancer_qc.h5ad"
QUICK_TEST_CELLS = 5_000

# GEO download URL for GSE176078 supplementary count matrix
GEO_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/"
    "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"
)

# ──────────────────────────────────────────────────────────────────────────────
# Breast cancer cell-type marker definitions (used for synthetic fallback &
# downstream annotation)
# ──────────────────────────────────────────────────────────────────────────────
CELL_TYPE_CONFIG = {
    "Cancer_Epithelial": {
        "fraction": 0.30,
        "markers": ["EPCAM", "KRT19", "KRT8", "KRT18", "MUC1", "ESR1",
                     "ERBB2", "MKI67", "TOP2A", "CD24"],
        "n_specific_genes": 200,
    },
    "T_cells": {
        "fraction": 0.18,
        "markers": ["CD3D", "CD3E", "CD3G", "CD2", "PTPRC", "IL7R",
                     "CD8A", "CD8B", "PDCD1", "HAVCR2"],
        "n_specific_genes": 120,
    },
    "B_cells": {
        "fraction": 0.08,
        "markers": ["CD79A", "CD79B", "MS4A1", "CD19", "BANK1", "PTPRC",
                     "PAX5", "IGHM", "IGHD", "TCL1A"],
        "n_specific_genes": 100,
    },
    "Macrophages": {
        "fraction": 0.12,
        "markers": ["CD68", "CD163", "CSF1R", "MRC1", "MSR1", "PTPRC",
                     "MARCO", "C1QA", "C1QB", "APOE"],
        "n_specific_genes": 120,
    },
    "Fibroblasts": {
        "fraction": 0.10,
        "markers": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP",
                     "ACTA2", "PDGFRA", "THY1", "TAGLN"],
        "n_specific_genes": 120,
    },
    "Endothelial": {
        "fraction": 0.06,
        "markers": ["PECAM1", "VWF", "CDH5", "CLDN5", "FLT1", "KDR",
                     "EMCN", "ERG", "ENG", "PLVAP"],
        "n_specific_genes": 100,
    },
    "NK_cells": {
        "fraction": 0.05,
        "markers": ["NKG7", "GNLY", "KLRD1", "KLRB1", "NCAM1", "PTPRC",
                     "PRF1", "GZMA", "GZMB", "FCGR3A"],
        "n_specific_genes": 80,
    },
    "Dendritic_cells": {
        "fraction": 0.04,
        "markers": ["FCER1A", "CD1C", "CLEC10A", "ITGAX", "HLA-DRA",
                     "HLA-DQA1", "PTPRC", "IRF8", "BATF3", "CLEC9A"],
        "n_specific_genes": 80,
    },
    "Mast_cells": {
        "fraction": 0.03,
        "markers": ["TPSAB1", "TPSB2", "KIT", "CPA3", "HPGDS", "PTPRC",
                     "HDC", "GATA2", "IL1RL1", "MS4A2"],
        "n_specific_genes": 60,
    },
    "Plasmablasts": {
        "fraction": 0.04,
        "markers": ["MZB1", "SDC1", "TNFRSF17", "XBP1", "PRDM1",
                     "PTPRC", "IRF4", "JCHAIN", "IGHA1", "IGHG1"],
        "n_specific_genes": 80,
    },
}


# ──────────────────────────────────────────────────────────────────────────────
# Utility helpers
# ──────────────────────────────────────────────────────────────────────────────

def _ensure_dirs() -> None:
    """Create required project directories."""
    for d in [RAW_DIR, PROCESSED_DIR, RESULTS_DIR]:
        d.mkdir(parents=True, exist_ok=True)


def _download_with_progress(url: str, dest: Path, desc: str = "Downloading") -> bool:
    """Download a file with a tqdm progress bar. Returns True on success."""
    try:
        import requests
        from tqdm import tqdm

        print(f"[INFO] {desc}: {url}")
        resp = requests.get(url, stream=True, timeout=60)
        resp.raise_for_status()
        total = int(resp.headers.get("content-length", 0))

        with open(dest, "wb") as f, tqdm(
            total=total, unit="B", unit_scale=True, desc=desc
        ) as pbar:
            for chunk in resp.iter_content(chunk_size=1 << 20):
                f.write(chunk)
                pbar.update(len(chunk))

        print(f"[INFO] Saved to {dest}  ({dest.stat().st_size / 1e6:.1f} MB)")
        return True

    except Exception as exc:
        print(f"[WARN] Download failed: {exc}")
        return False


# ──────────────────────────────────────────────────────────────────────────────
# Strategy 1: Download from GEO
# ──────────────────────────────────────────────────────────────────────────────

def _load_from_geo() -> ad.AnnData | None:
    """Try to download and parse the GSE176078 count matrix from GEO/NCBI."""
    archive = RAW_DIR / "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"

    # Download if not cached
    if not archive.exists():
        ok = _download_with_progress(GEO_URL, archive, "Downloading GSE176078")
        if not ok:
            return None

    print("[INFO] Extracting count matrices from archive …")
    try:
        with tarfile.open(archive, "r:gz") as tar:
            if sys.version_info >= (3, 12):
                tar.extractall(RAW_DIR, filter="data")
            else:
                tar.extractall(RAW_DIR)
    except Exception as exc:
        print(f"[WARN] Extraction failed: {exc}")
        return None

    # Look for 10x-format directories (barcodes.tsv, features/genes.tsv, matrix.mtx)
    adatas = []
    for mtx_file in sorted(RAW_DIR.rglob("matrix.mtx*")):
        sample_dir = mtx_file.parent
        sample_name = sample_dir.name
        try:
            adata_sample = sc.read_10x_mtx(sample_dir, var_names="gene_symbols")
            adata_sample.obs["sample"] = sample_name
            adata_sample.var_names_make_unique()
            adatas.append(adata_sample)
            print(f"  ✓ Loaded {sample_name}: {adata_sample.n_obs} cells")
        except Exception as e:
            print(f"  ✗ Skipped {sample_name}: {e}")

    if not adatas:
        # Try GSE176078-style flat layout: count_matrix_sparse.mtx + genes/barcodes TSVs
        mtx_file = next(RAW_DIR.rglob("count_matrix_sparse.mtx*"), None)
        genes_file = next(RAW_DIR.rglob("count_matrix_genes.tsv*"), None)
        barcodes_file = next(RAW_DIR.rglob("count_matrix_barcodes.tsv*"), None)
        if mtx_file and genes_file and barcodes_file:
            try:
                from scipy.io import mmread
                mat = mmread(str(mtx_file)).T  # MTX is genes×cells, transpose to cells×genes
                genes = pd.read_csv(genes_file, header=None, sep="\t")[0].values
                barcodes = pd.read_csv(barcodes_file, header=None, sep="\t")[0].values
                adata = ad.AnnData(
                    sp.csr_matrix(mat),
                    obs=pd.DataFrame(index=barcodes),
                    var=pd.DataFrame(index=genes),
                )
                adata.var_names_make_unique()
                adatas.append(adata)
                print(f"  ✓ Loaded sparse matrix: {adata.n_obs} cells × {adata.n_vars} genes")
            except Exception as e:
                print(f"  ✗ Failed to load sparse matrix: {e}")

    if not adatas:
        # Try reading as a single CSV/TSV if 10x format not found
        for csv_file in sorted(RAW_DIR.rglob("*.csv*")) + sorted(RAW_DIR.rglob("*.tsv*")):
            try:
                sep = "\t" if "tsv" in csv_file.name else ","
                df = pd.read_csv(csv_file, index_col=0, sep=sep, nrows=5)
                print(f"  Found tabular file: {csv_file.name}")
                df_full = pd.read_csv(csv_file, index_col=0, sep=sep)
                adata = ad.AnnData(sp.csr_matrix(df_full.values),
                                   obs=pd.DataFrame(index=df_full.index),
                                   var=pd.DataFrame(index=df_full.columns))
                adatas.append(adata)
                break
            except Exception:
                continue

    if not adatas:
        print("[WARN] Could not parse any count matrices from archive.")
        return None

    # Concatenate all samples
    if len(adatas) == 1:
        adata = adatas[0]
    else:
        adata = ad.concat(adatas, join="outer", label="sample", index_unique="-")

    adata.var_names_make_unique()
    print(f"[INFO] Loaded GEO data: {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


# ──────────────────────────────────────────────────────────────────────────────
# Strategy 2: Generate high-fidelity synthetic dataset
# ──────────────────────────────────────────────────────────────────────────────

def _generate_synthetic_breast_cancer(
    n_cells: int = 20_000,
    n_genes: int = 5_000,
    seed: int = 42,
) -> ad.AnnData:
    """
    Generate a realistic synthetic breast cancer scRNA-seq dataset.

    Cell-type proportions, marker gene expression patterns, and count
    distributions are modelled after published breast cancer atlases
    (Wu et al. 2021, Pal et al. 2021).  Counts follow a negative-binomial
    distribution with cell-type-specific means and dispersions.
    """
    rng = np.random.default_rng(seed)
    print(f"[INFO] Generating synthetic breast cancer dataset "
          f"({n_cells} cells, {n_genes} genes) …")

    # ------------------------------------------------------------------
    # 1. Assign cells to cell types
    # ------------------------------------------------------------------
    cell_types = list(CELL_TYPE_CONFIG.keys())
    fractions = np.array([CELL_TYPE_CONFIG[ct]["fraction"] for ct in cell_types])
    fractions /= fractions.sum()  # normalize
    labels = rng.choice(cell_types, size=n_cells, p=fractions)

    # ------------------------------------------------------------------
    # 2. Build gene names: markers + cell-type-specific + housekeeping
    # ------------------------------------------------------------------
    all_markers: list[str] = []
    for cfg in CELL_TYPE_CONFIG.values():
        all_markers.extend(cfg["markers"])
    all_markers = list(dict.fromkeys(all_markers))  # deduplicate, keep order

    # Remaining genes: generic gene names
    n_remaining = n_genes - len(all_markers)
    generic_genes = [f"GENE_{i:04d}" for i in range(n_remaining)]
    gene_names = all_markers + generic_genes

    # ------------------------------------------------------------------
    # 3. Generate count matrix (negative-binomial)
    # ------------------------------------------------------------------
    # Base expression: low background for all genes
    base_mean = 0.3
    base_disp = 0.5  # dispersion (smaller = more variable)

    # Pre-compute prob/n for NB parameterization: mean = n*p/(1-p)
    def _nb_params(mean: float, disp: float):
        """Convert mean+dispersion to scipy NB parameters."""
        p = disp / (disp + mean)
        return disp, p

    n_nb, p_nb = _nb_params(base_mean, base_disp)
    X = np.zeros((n_cells, n_genes), dtype=np.float32)

    # Background expression
    X[:] = rng.negative_binomial(max(1, int(n_nb)), min(p_nb, 0.999),
                                  size=(n_cells, n_genes)).astype(np.float32)

    # ------------------------------------------------------------------
    # 4. Inject cell-type-specific marker expression
    # ------------------------------------------------------------------
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    for ct_name, cfg in CELL_TYPE_CONFIG.items():
        ct_mask = labels == ct_name
        if ct_mask.sum() == 0:
            continue

        # Strong marker expression in correct cell type
        for marker in cfg["markers"]:
            if marker in gene_to_idx:
                idx = gene_to_idx[marker]
                marker_mean = rng.uniform(4.0, 12.0)  # high expression
                n_m, p_m = _nb_params(marker_mean, 2.0)
                X[ct_mask, idx] = rng.negative_binomial(
                    max(1, int(n_m)), min(p_m, 0.999), size=ct_mask.sum()
                )

        # Moderate expression for cell-type-specific generic genes
        specific_start = len(all_markers)
        ct_idx = cell_types.index(ct_name)
        block_size = cfg["n_specific_genes"]
        block_start = specific_start + ct_idx * 200  # non-overlapping blocks
        block_end = min(block_start + block_size, n_genes)

        if block_start < n_genes:
            specific_mean = rng.uniform(2.0, 5.0)
            n_s, p_s = _nb_params(specific_mean, 1.5)
            X[ct_mask, block_start:block_end] = rng.negative_binomial(
                max(1, int(n_s)), min(p_s, 0.999),
                size=(ct_mask.sum(), block_end - block_start),
            )

    # ------------------------------------------------------------------
    # 5. Add housekeeping genes (expressed in all cells)
    # ------------------------------------------------------------------
    housekeeping = ["ACTB", "GAPDH", "B2M", "MALAT1", "TMSB4X"]
    for hk in housekeeping:
        if hk not in gene_to_idx:
            # Replace a generic gene
            replace_idx = len(all_markers) + len(generic_genes) - 1
            if replace_idx < n_genes:
                gene_names[replace_idx] = hk
                gene_to_idx[hk] = replace_idx
                generic_genes.pop()

        if hk in gene_to_idx:
            idx = gene_to_idx[hk]
            hk_mean = rng.uniform(5.0, 15.0)
            n_h, p_h = _nb_params(hk_mean, 3.0)
            X[:, idx] = rng.negative_binomial(
                max(1, int(n_h)), min(p_h, 0.999), size=n_cells
            )

    # ------------------------------------------------------------------
    # 6. Add mitochondrial genes (for QC)
    # ------------------------------------------------------------------
    mito_genes = ["MT-CO1", "MT-CO2", "MT-CO3", "MT-ND1", "MT-ND4",
                  "MT-ATP6", "MT-CYB", "MT-ND5", "MT-ND2", "MT-RNR2"]
    for i, mg in enumerate(mito_genes):
        replace_idx = len(all_markers) + i
        if replace_idx < n_genes:
            gene_names[replace_idx] = mg
            gene_to_idx[mg] = replace_idx
            # Low-moderate mito expression
            mito_mean = rng.uniform(1.0, 4.0)
            n_mt, p_mt = _nb_params(mito_mean, 1.0)
            X[:, replace_idx] = rng.negative_binomial(
                max(1, int(n_mt)), min(p_mt, 0.999), size=n_cells
            )

    # Spike a few cells with high mito (to be filtered in QC)
    n_high_mito = int(0.03 * n_cells)
    high_mito_cells = rng.choice(n_cells, n_high_mito, replace=False)
    for mg in mito_genes:
        if mg in gene_to_idx:
            idx = gene_to_idx[mg]
            X[high_mito_cells, idx] = rng.negative_binomial(20, 0.3,
                                                             size=n_high_mito)

    # ------------------------------------------------------------------
    # 7. Simulate sample / patient metadata
    # ------------------------------------------------------------------
    n_patients = 5
    patient_ids = [f"Patient_{i+1}" for i in range(n_patients)]
    conditions = ["Tumor", "Tumor", "Tumor", "Normal_adjacent", "Normal_adjacent"]
    subtypes = ["TNBC", "Luminal_A", "HER2+", "Normal", "Normal"]

    patient_labels = rng.choice(patient_ids, size=n_cells)
    condition_labels = [conditions[patient_ids.index(p)] for p in patient_labels]
    subtype_labels = [subtypes[patient_ids.index(p)] for p in patient_labels]

    # ------------------------------------------------------------------
    # 8. Assemble AnnData object
    # ------------------------------------------------------------------
    adata = ad.AnnData(
        X=sp.csr_matrix(X),
        obs=pd.DataFrame(
            {
                "cell_type_ground_truth": pd.Categorical(labels),
                "patient_id": pd.Categorical(patient_labels),
                "condition": pd.Categorical(condition_labels),
                "subtype": pd.Categorical(subtype_labels),
            },
            index=[f"Cell_{i:06d}" for i in range(n_cells)],
        ),
        var=pd.DataFrame(index=gene_names),
    )
    adata.var_names_make_unique()

    # Add some batch effect per patient (multiplicative scaling)
    for pid in patient_ids:
        mask = (adata.obs["patient_id"] == pid).values
        scale = rng.uniform(0.8, 1.3)
        adata.X[mask] = (adata.X[mask].toarray() * scale).astype(np.float32)
        adata.X[mask] = sp.csr_matrix(adata.X[mask])

    adata.X = sp.csr_matrix(np.round(adata.X.toarray()).clip(0).astype(np.float32))

    print(f"[INFO] Synthetic dataset ready: {adata.n_obs} cells × {adata.n_vars} genes")
    print(f"       Cell types: {dict(zip(*np.unique(labels, return_counts=True)))}")
    return adata


# ──────────────────────────────────────────────────────────────────────────────
# Quality control pipeline
# ──────────────────────────────────────────────────────────────────────────────

def run_qc(adata: ad.AnnData) -> ad.AnnData:
    """Run standard scRNA-seq quality control filters and metrics."""
    print("\n" + "=" * 60)
    print("QUALITY CONTROL")
    print("=" * 60)

    n_before = adata.n_obs
    print(f"\n[QC] Starting cells: {adata.n_obs:,}")
    print(f"[QC] Starting genes: {adata.n_vars:,}")

    # ── Compute QC metrics ────────────────────────────────────────────────
    # Mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # Ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # Hemoglobin genes
    adata.var["hb"] = adata.var_names.str.match("^HB[^(P)]")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], percent_top=None,
        log1p=False, inplace=True,
    )

    print(f"\n[QC] Median genes/cell:  {adata.obs['n_genes_by_counts'].median():.0f}")
    print(f"[QC] Median UMIs/cell:   {adata.obs['total_counts'].median():.0f}")
    print(f"[QC] Median %mito:       {adata.obs['pct_counts_mt'].median():.2f}%")

    # ── QC violin plots ──────────────────────────────────────────────────
    try:
        fig, axes = __import__("matplotlib.pyplot", fromlist=["pyplot"]).subplots(
            1, 4, figsize=(16, 4)
        )
        sc.pl.violin(adata, "n_genes_by_counts", ax=axes[0], show=False)
        sc.pl.violin(adata, "total_counts", ax=axes[1], show=False)
        sc.pl.violin(adata, "pct_counts_mt", ax=axes[2], show=False)
        sc.pl.violin(adata, "pct_counts_ribo", ax=axes[3], show=False)
        fig.tight_layout()
        fig.savefig(RESULTS_DIR / "qc_violins.png", dpi=150, bbox_inches="tight")
        __import__("matplotlib.pyplot", fromlist=["pyplot"]).close(fig)
        print("[QC] Saved QC violin plot → results/qc_violins.png")
    except Exception as exc:
        print(f"[QC] Could not save QC plot: {exc}")

    # ── Filter cells ─────────────────────────────────────────────────────
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_cells(adata, max_genes=6000)
    adata = adata[adata.obs["pct_counts_mt"] < 20, :].copy()

    # ── Filter genes ─────────────────────────────────────────────────────
    sc.pp.filter_genes(adata, min_cells=10)

    n_after = adata.n_obs
    print(f"\n[QC] Cells after filtering: {n_after:,}  "
          f"(removed {n_before - n_after:,}, {(n_before-n_after)/n_before*100:.1f}%)")
    print(f"[QC] Genes after filtering: {adata.n_vars:,}")

    return adata


# ──────────────────────────────────────────────────────────────────────────────
# Main entry point
# ──────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download breast cancer scRNA-seq data and run QC."
    )
    parser.add_argument(
        "--mode",
        choices=["full", "quick-test"],
        default="quick-test",
        help="'full' for all cells; 'quick-test' for ~5,000 cells (default).",
    )
    parser.add_argument(
        "--force-synthetic",
        action="store_true",
        help="Skip download and generate synthetic data directly.",
    )
    args = parser.parse_args()

    _ensure_dirs()

    # Check for cached data
    if QC_H5AD.exists():
        print(f"[INFO] Found cached QC data: {QC_H5AD}")
        adata = sc.read_h5ad(QC_H5AD)
        print(f"[INFO] Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
        print("[INFO] Delete the cached file to re-run QC from scratch.")
        return

    # ── Load data ─────────────────────────────────────────────────────────
    adata = None

    if not args.force_synthetic:
        print("\n[STEP 1] Attempting to download from GEO (GSE176078) …")
        adata = _load_from_geo()

    if adata is None:
        print("\n[STEP 1] Using synthetic breast cancer dataset (fallback) …")
        n_cells = 20_000 if args.mode == "full" else QUICK_TEST_CELLS + 1_000
        adata = _generate_synthetic_breast_cancer(n_cells=n_cells)

    # ── Subset for quick-test mode ────────────────────────────────────────
    if args.mode == "quick-test" and adata.n_obs > QUICK_TEST_CELLS:
        print(f"\n[INFO] Quick-test mode: subsetting to {QUICK_TEST_CELLS} cells …")
        sc.pp.subsample(adata, n_obs=QUICK_TEST_CELLS, random_state=42)
        print(f"[INFO] Subsetted to {adata.n_obs} cells")

    # ── Run QC ────────────────────────────────────────────────────────────
    adata = run_qc(adata)

    # ── Store raw counts before normalization (for DE analysis later) ─────
    adata.layers["raw_counts"] = adata.X.copy()

    # ── Save ──────────────────────────────────────────────────────────────
    adata.write_h5ad(QC_H5AD)
    print(f"\n[INFO] Saved QC-filtered data → {QC_H5AD}")
    print(f"[INFO] Final: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    print("\n✅ Step 1 complete. Run 02_clustering_and_annotation.py next.")


if __name__ == "__main__":
    main()
