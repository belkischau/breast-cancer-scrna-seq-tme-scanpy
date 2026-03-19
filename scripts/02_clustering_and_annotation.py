#!/usr/bin/env python3
"""
02_clustering_and_annotation.py
===============================
Normalize, find HVGs, reduce dimensions, cluster, and annotate cell types.

Reads:   data/processed/breast_cancer_qc.h5ad
Writes:  data/processed/breast_cancer_annotated.h5ad
         results/umap_celltype.png
         results/umap_patient.png
         results/umap_condition.png

Usage:
    python scripts/02_clustering_and_annotation.py
    python scripts/02_clustering_and_annotation.py --resolution 1.0
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

warnings.filterwarnings("ignore", category=FutureWarning)

# ──────────────────────────────────────────────────────────────────────────────
# Paths
# ──────────────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"

INPUT_H5AD = PROCESSED_DIR / "breast_cancer_qc.h5ad"
OUTPUT_H5AD = PROCESSED_DIR / "breast_cancer_annotated.h5ad"

# ──────────────────────────────────────────────────────────────────────────────
# Canonical marker genes for breast cancer TME cell types
# ──────────────────────────────────────────────────────────────────────────────
MARKER_GENES = {
    "Cancer_Epithelial": ["EPCAM", "KRT19", "KRT8", "KRT18", "MUC1", "CD24"],
    "T_cells":           ["CD3D", "CD3E", "CD2", "IL7R", "CD8A"],
    "B_cells":           ["CD79A", "MS4A1", "CD19", "BANK1"],
    "Macrophages":       ["CD68", "CD163", "CSF1R", "MRC1", "C1QA"],
    "Fibroblasts":       ["COL1A1", "COL1A2", "COL3A1", "DCN", "FAP"],
    "Endothelial":       ["PECAM1", "VWF", "CDH5", "CLDN5"],
    "NK_cells":          ["NKG7", "GNLY", "KLRD1", "PRF1"],
    "Dendritic_cells":   ["FCER1A", "CD1C", "CLEC10A", "HLA-DRA"],
    "Mast_cells":        ["TPSAB1", "KIT", "CPA3", "HPGDS"],
    "Plasmablasts":      ["MZB1", "SDC1", "XBP1", "JCHAIN"],
}

# Immune marker for broad annotation
IMMUNE_MARKER = "PTPRC"  # CD45


# ──────────────────────────────────────────────────────────────────────────────
# Preprocessing & dimensionality reduction
# ──────────────────────────────────────────────────────────────────────────────

def preprocess(adata: ad.AnnData) -> ad.AnnData:
    """Normalize, log-transform, select HVGs, scale, PCA, neighbors, UMAP."""
    print("\n" + "=" * 60)
    print("PREPROCESSING & DIMENSIONALITY REDUCTION")
    print("=" * 60)

    # ── Normalization ─────────────────────────────────────────────────────
    print("\n[STEP] Library-size normalization (10,000 CPM) + log1p …")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Store normalized counts before HVG + scaling (for visualization later)
    adata.raw = adata.copy()

    # ── Highly Variable Genes ─────────────────────────────────────────────
    print("[STEP] Selecting top 2,000 highly variable genes (Pearson residuals) …")
    sc.experimental.pp.highly_variable_genes(
        adata, n_top_genes=2000, flavor="pearson_residuals",
        subset=False, layer="raw_counts",
    )
    n_hvg = adata.var["highly_variable"].sum()
    print(f"       Found {n_hvg} HVGs")

    # Subset to HVGs for downstream analysis
    adata = adata[:, adata.var["highly_variable"]].copy()

    # ── Scaling ───────────────────────────────────────────────────────────
    print("[STEP] Scaling to unit variance (max_value=10) …")
    sc.pp.scale(adata, max_value=10)

    # ── PCA ───────────────────────────────────────────────────────────────
    n_pcs = min(50, adata.n_obs - 1, adata.n_vars - 1)
    print(f"[STEP] PCA ({n_pcs} components) …")
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")

    # Variance explained
    cumvar = np.cumsum(adata.uns["pca"]["variance_ratio"])
    n_90 = int(np.searchsorted(cumvar, 0.90)) + 1
    print(f"       {n_90} PCs explain 90% of variance")

    # ── Neighbors + UMAP ─────────────────────────────────────────────────
    n_neighbors = min(30, adata.n_obs - 1)
    print(f"[STEP] Computing neighbors (k={n_neighbors}) and UMAP …")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=min(30, n_pcs))
    sc.tl.umap(adata, random_state=42)

    return adata


# ──────────────────────────────────────────────────────────────────────────────
# Clustering
# ──────────────────────────────────────────────────────────────────────────────

def cluster(adata: ad.AnnData, resolution: float = 1.0) -> ad.AnnData:
    """Run Leiden clustering."""
    print(f"\n[STEP] Leiden clustering (resolution={resolution}) …")
    sc.tl.leiden(adata, resolution=resolution, random_state=42, flavor="igraph",
                 n_iterations=2)

    n_clusters = adata.obs["leiden"].nunique()
    print(f"       Found {n_clusters} clusters")

    # Cluster sizes
    cluster_counts = adata.obs["leiden"].value_counts().sort_index()
    for cl, cnt in cluster_counts.items():
        print(f"       Cluster {cl}: {cnt:,} cells ({cnt/adata.n_obs*100:.1f}%)")

    return adata


# ──────────────────────────────────────────────────────────────────────────────
# Cell-type annotation
# ──────────────────────────────────────────────────────────────────────────────

def annotate_cell_types(adata: ad.AnnData) -> ad.AnnData:
    """
    Annotate clusters to cell types using marker gene scoring.

    For each cluster, compute the mean expression of each cell-type's marker
    genes (from .raw) and assign the cell type with the highest score.
    """
    print("\n" + "=" * 60)
    print("CELL TYPE ANNOTATION")
    print("=" * 60)

    # Use raw (normalized, log-transformed) expression for scoring
    if adata.raw is not None:
        raw_adata = adata.raw.to_adata()
    else:
        raw_adata = adata

    # Build a score matrix: clusters × cell_types
    clusters = sorted(adata.obs["leiden"].unique(), key=int)
    scores = pd.DataFrame(0.0, index=clusters, columns=list(MARKER_GENES.keys()))

    for ct_name, markers in MARKER_GENES.items():
        # Keep only markers present in the dataset
        available = [m for m in markers if m in raw_adata.var_names]
        if not available:
            print(f"  [WARN] No markers found for {ct_name}")
            continue

        # Score each cell
        sc.tl.score_genes(adata, gene_list=available, score_name=f"score_{ct_name}",
                          use_raw=True)

        # Aggregate score per cluster
        for cl in clusters:
            mask = adata.obs["leiden"] == cl
            scores.loc[cl, ct_name] = adata.obs.loc[mask, f"score_{ct_name}"].mean()

    # Assign best-matching cell type to each cluster
    cluster_to_ct = scores.idxmax(axis=1).to_dict()

    # Handle ties / duplicates: append cluster number if same type assigned twice
    ct_counts: dict[str, int] = {}
    final_mapping: dict[str, str] = {}
    for cl in clusters:
        ct = cluster_to_ct[cl]
        ct_counts[ct] = ct_counts.get(ct, 0) + 1
        if ct_counts[ct] > 1:
            final_mapping[cl] = f"{ct}_{ct_counts[ct]}"
        else:
            final_mapping[cl] = ct

    adata.obs["cell_type"] = (
        adata.obs["leiden"].map(final_mapping).astype("category")
    )

    print("\n  Cluster → Cell Type mapping:")
    print("  " + "-" * 40)
    for cl, ct in sorted(final_mapping.items(), key=lambda x: int(x[0])):
        n = (adata.obs["leiden"] == cl).sum()
        top_score = scores.loc[cl, cluster_to_ct[cl]]
        print(f"  Cluster {cl:>2s} → {ct:<25s} ({n:,} cells, score={top_score:.3f})")

    # Summary
    print(f"\n  Total cell types: {adata.obs['cell_type'].nunique()}")
    print("\n  Cell type distribution:")
    for ct, cnt in adata.obs["cell_type"].value_counts().items():
        print(f"    {ct}: {cnt:,} ({cnt/adata.n_obs*100:.1f}%)")

    return adata


# ──────────────────────────────────────────────────────────────────────────────
# Marker verification
# ──────────────────────────────────────────────────────────────────────────────

def verify_markers(adata: ad.AnnData) -> None:
    """Verify that canonical markers are expressed in the expected cell types."""
    print("\n" + "=" * 60)
    print("MARKER VERIFICATION")
    print("=" * 60)

    key_checks = {
        "EPCAM":   "Cancer_Epithelial",
        "KRT19":   "Cancer_Epithelial",
        "PTPRC":   "T_cells",       # CD45 — should be high in all immune
        "CD3D":    "T_cells",
        "CD79A":   "B_cells",
        "CD68":    "Macrophages",
        "COL1A1":  "Fibroblasts",
        "PECAM1":  "Endothelial",
        "NKG7":    "NK_cells",
        "TPSAB1":  "Mast_cells",
    }

    raw_adata = adata.raw.to_adata() if adata.raw is not None else adata
    passed = 0
    total = 0

    for gene, expected_ct in key_checks.items():
        total += 1
        if gene not in raw_adata.var_names:
            print(f"  ⚠  {gene:<10s} — not in dataset, skipping")
            continue

        # Mean expression per cell type
        expr = pd.Series(
            raw_adata[:, gene].X.toarray().ravel()
            if hasattr(raw_adata[:, gene].X, "toarray")
            else raw_adata[:, gene].X.ravel(),
            index=adata.obs.index,
        )
        mean_per_ct = expr.groupby(adata.obs["cell_type"]).mean()

        top_ct = mean_per_ct.idxmax()
        # Check if the expected cell type is the top or among the top 2
        top2 = mean_per_ct.nlargest(2).index.tolist()
        # Also match partial names (e.g. "Cancer_Epithelial_2" matches "Cancer_Epithelial")
        match = any(expected_ct in ct for ct in top2)

        if match:
            print(f"  ✅ {gene:<10s} highest in {top_ct:<25s} (expected: {expected_ct})")
            passed += 1
        else:
            print(f"  ❌ {gene:<10s} highest in {top_ct:<25s} (expected: {expected_ct})")

    print(f"\n  Verification: {passed}/{total} markers in expected cell types")


# ──────────────────────────────────────────────────────────────────────────────
# Plotting
# ──────────────────────────────────────────────────────────────────────────────

def plot_umaps(adata: ad.AnnData) -> None:
    """Generate UMAP plots colored by cell type, patient, and condition."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    sc.set_figure_params(dpi=150, frameon=False, figsize=(8, 6))

    # ── UMAP by cell type ─────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color="cell_type", ax=ax, show=False, legend_loc="on data",
               legend_fontsize=7, title="Breast Cancer TME — Cell Types")
    fig.savefig(RESULTS_DIR / "umap_celltype.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[PLOT] Saved → results/umap_celltype.png")

    # ── UMAP by Leiden cluster ────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color="leiden", ax=ax, show=False, legend_loc="on data",
               title="Leiden Clusters")
    fig.savefig(RESULTS_DIR / "umap_leiden.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[PLOT] Saved → results/umap_leiden.png")

    # ── UMAP by patient ──────────────────────────────────────────────────
    if "patient_id" not in adata.obs.columns:
        # Extract patient ID from barcode prefix (e.g., "CID3586_AAGACCTCAGCATGAG" → "CID3586")
        if adata.obs.index.str.contains("_").any():
            adata.obs["patient_id"] = adata.obs.index.str.split("_").str[0]

    if "patient_id" in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(10, 8))
        sc.pl.umap(adata, color="patient_id", ax=ax, show=False,
                    title="Patient ID")
        fig.savefig(RESULTS_DIR / "umap_patient.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"[PLOT] Saved → results/umap_patient.png")

    # ── UMAP by condition ────────────────────────────────────────────────
    if "condition" in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(10, 8))
        sc.pl.umap(adata, color="condition", ax=ax, show=False,
                    title="Condition (Tumor vs. Normal)")
        fig.savefig(RESULTS_DIR / "umap_condition.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"[PLOT] Saved → results/umap_condition.png")


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="Cluster and annotate cell types.")
    parser.add_argument("--resolution", type=float, default=1.0,
                        help="Leiden clustering resolution (default: 1.0).")
    args = parser.parse_args()

    if not INPUT_H5AD.exists():
        print(f"[ERROR] Input not found: {INPUT_H5AD}")
        print("        Run 01_download_and_qc.py first.")
        return

    # Load
    print(f"[INFO] Loading {INPUT_H5AD} …")
    adata = sc.read_h5ad(INPUT_H5AD)
    print(f"[INFO] Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    # Pipeline
    adata = preprocess(adata)
    adata = cluster(adata, resolution=args.resolution)
    adata = annotate_cell_types(adata)
    verify_markers(adata)
    plot_umaps(adata)

    # Save
    adata.write_h5ad(OUTPUT_H5AD)
    print(f"\n[INFO] Saved annotated data → {OUTPUT_H5AD}")
    print("\n✅ Step 2 complete. Run 03_de_analysis_and_markers.py next.")


if __name__ == "__main__":
    main()
