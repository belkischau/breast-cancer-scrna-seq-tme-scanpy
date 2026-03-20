#!/usr/bin/env python3
"""
03_de_analysis_and_markers.py
=============================
Differential expression analysis and marker gene identification.

Analyses performed:
  1. Marker genes per cell type (one-vs-rest Wilcoxon rank-sum)
  2. Cancer-specific DEGs (cancer epithelial vs. all other cell types)
  3. Condition-specific DEGs (tumor vs. normal within cancer epithelial)
  4. Volcano plots for top DEGs
  5. Marker validation with biological context

Reads:   data/processed/breast_cancer_annotated.h5ad
Writes:  data/processed/breast_cancer_de.h5ad
         results/deg_volcano.png
         results/marker_dotplot.png
         results/top_markers_per_celltype.csv

Usage:
    python scripts/03_de_analysis_and_markers.py
"""

from __future__ import annotations

import warnings
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

# ──────────────────────────────────────────────────────────────────────────────
# Paths
# ──────────────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"

INPUT_H5AD = PROCESSED_DIR / "breast_cancer_annotated.h5ad"
OUTPUT_H5AD = PROCESSED_DIR / "breast_cancer_de.h5ad"

# ──────────────────────────────────────────────────────────────────────────────
# Biological context for key DEGs
# ──────────────────────────────────────────────────────────────────────────────
GENE_BIO_CONTEXT = {
    "EPCAM": "Epithelial cell adhesion molecule — cancer epithelial marker",
    "KRT19": "Cytokeratin 19 — luminal breast cancer marker",
    "KRT8": "Cytokeratin 8 — epithelial differentiation",
    "MUC1": "Mucin-1 — overexpressed in breast cancer, immunotherapy target",
    "MKI67": "Ki-67 — proliferation marker, prognostic indicator",
    "TOP2A": "Topoisomerase IIα — drug target (anthracyclines)",
    "ESR1": "Estrogen receptor α — hormone therapy target",
    "ERBB2": "HER2 — amplified in HER2+ breast cancer, trastuzumab target",
    "CD274": "PD-L1 — immune checkpoint, immunotherapy biomarker",
    "CD3D": "T cell receptor complex — pan-T cell marker",
    "CD8A": "Cytotoxic T cell marker",
    "PDCD1": "PD-1 — T cell exhaustion / immune checkpoint",
    "HAVCR2": "TIM-3 — T cell exhaustion marker",
    "LAG3": "LAG-3 — immune checkpoint receptor",
    "CD68": "Macrophage marker",
    "CD163": "M2 macrophage marker — immunosuppressive TAMs",
    "COL1A1": "Collagen I — CAF / stromal marker",
    "FAP": "Fibroblast activation protein — CAF marker",
    "PECAM1": "CD31 — endothelial marker, angiogenesis",
    "VEGFA": "VEGF-A — angiogenesis driver",
    "NKG7": "NK cell granule protein — cytotoxicity",
    "FOXP3": "Regulatory T cell master transcription factor",
    "IL6": "Inflammatory cytokine — iCAF marker",
    "ACTA2": "α-SMA — myofibroblast / myCAF marker",
}


# ──────────────────────────────────────────────────────────────────────────────
# 1. Marker genes per cell type
# ──────────────────────────────────────────────────────────────────────────────

def find_cell_type_markers(adata: ad.AnnData) -> pd.DataFrame:
    """
    Run one-vs-rest Wilcoxon rank-sum test to find marker genes for each
    annotated cell type.
    """
    print("\n" + "=" * 60)
    print("MARKER GENE IDENTIFICATION (one-vs-rest)")
    print("=" * 60)

    sc.tl.rank_genes_groups(
        adata,
        groupby="cell_type",
        method="wilcoxon",
        use_raw=True,
        pts=True,
        key_added="markers_celltype",
    )

    # Extract results into a tidy DataFrame
    result = sc.get.rank_genes_groups_df(adata, group=None, key="markers_celltype")
    print(f"\n  Total marker tests: {len(result):,}")

    # Filter significant
    sig = result[result["pvals_adj"] < 0.05].copy()
    print(f"  Significant (FDR < 0.05): {sig.shape[0]:,}")

    # Top 10 per cell type
    print("\n  ── Top 10 markers per cell type ─────────────────────────")
    for ct in sorted(adata.obs["cell_type"].unique()):
        ct_markers = sig[sig["group"] == ct].nlargest(10, "scores")
        if ct_markers.empty:
            continue
        genes = ct_markers["names"].tolist()
        scores = ct_markers["scores"].tolist()
        print(f"\n  {ct}:")
        for g, s in zip(genes, scores):
            bio = GENE_BIO_CONTEXT.get(g, "")
            bio_str = f"  — {bio}" if bio else ""
            print(f"    {g:<12s}  score={s:>7.2f}{bio_str}")

    # Save to CSV
    top_markers = (
        sig.sort_values(["group", "scores"], ascending=[True, False])
        .groupby("group").head(20)
    )
    csv_path = RESULTS_DIR / "top_markers_per_celltype.csv"
    top_markers.to_csv(csv_path, index=False)
    print(f"\n  Saved top markers → {csv_path}")

    return sig


# ──────────────────────────────────────────────────────────────────────────────
# 2. Cancer-specific differential expression
# ──────────────────────────────────────────────────────────────────────────────

def cancer_vs_other_de(adata: ad.AnnData) -> pd.DataFrame | None:
    """
    Differential expression: Cancer epithelial cells vs. all other cell types.
    """
    print("\n" + "=" * 60)
    print("CANCER vs. OTHER CELL TYPES — DE ANALYSIS")
    print("=" * 60)

    # Identify cancer clusters
    cancer_types = [ct for ct in adata.obs["cell_type"].unique()
                    if "cancer" in ct.lower() or "epithelial" in ct.lower()]

    if not cancer_types:
        print("  [WARN] No cancer epithelial clusters found. Skipping.")
        return None

    print(f"  Cancer groups: {cancer_types}")

    # Create binary label
    adata.obs["is_cancer"] = (
        adata.obs["cell_type"].isin(cancer_types).map({True: "Cancer", False: "Other"})
    ).astype("category")

    sc.tl.rank_genes_groups(
        adata,
        groupby="is_cancer",
        groups=["Cancer"],
        reference="Other",
        method="wilcoxon",
        use_raw=True,
        key_added="cancer_vs_other",
    )

    result = sc.get.rank_genes_groups_df(adata, group="Cancer",
                                          key="cancer_vs_other")

    # Compute log2 fold change
    result["log2fc"] = result["logfoldchanges"]  # scanpy already gives log2FC

    sig_up = result[(result["pvals_adj"] < 0.05) & (result["logfoldchanges"] > 1)]
    sig_down = result[(result["pvals_adj"] < 0.05) & (result["logfoldchanges"] < -1)]

    print(f"\n  Total DEGs (FDR < 0.05): {len(result[result['pvals_adj'] < 0.05]):,}")
    print(f"  Upregulated in cancer (log2FC > 1): {len(sig_up):,}")
    print(f"  Downregulated in cancer (log2FC < -1): {len(sig_down):,}")

    # Top upregulated
    print("\n  ── Top 15 upregulated in cancer epithelial ──────────────")
    for _, row in sig_up.nlargest(15, "scores").iterrows():
        gene = row["names"]
        bio = GENE_BIO_CONTEXT.get(gene, "")
        bio_str = f"  — {bio}" if bio else ""
        print(f"    {gene:<12s} log2FC={row['logfoldchanges']:>6.2f}  "
              f"padj={row['pvals_adj']:.2e}{bio_str}")

    # Save
    csv_path = RESULTS_DIR / "cancer_vs_other_degs.csv"
    result.to_csv(csv_path, index=False)
    print(f"\n  Saved → {csv_path}")

    return result


# ──────────────────────────────────────────────────────────────────────────────
# 3. Condition-specific DE (Tumor vs Normal within cancer cells)
# ──────────────────────────────────────────────────────────────────────────────

def tumor_vs_normal_de(adata: ad.AnnData) -> pd.DataFrame | None:
    """DE between Tumor and Normal conditions within cancer epithelial cells."""
    print("\n" + "=" * 60)
    print("TUMOR vs. NORMAL (within cancer epithelial) — DE ANALYSIS")
    print("=" * 60)

    if "condition" not in adata.obs.columns:
        print("  [WARN] No 'condition' column found. Skipping.")
        return None

    # Subset to cancer epithelial
    cancer_mask = adata.obs["cell_type"].str.contains("Cancer|Epithelial",
                                                       case=False, na=False)
    adata_cancer = adata[cancer_mask].copy()

    conditions = adata_cancer.obs["condition"].unique()
    if len(conditions) < 2:
        print(f"  [WARN] Only one condition found ({conditions}). Skipping.")
        return None

    # Ensure we have Tumor and Normal-like groups
    tumor_group = [c for c in conditions if "tumor" in c.lower()]
    normal_group = [c for c in conditions if "normal" in c.lower()]

    if not tumor_group or not normal_group:
        print(f"  [WARN] Cannot identify Tumor/Normal groups in {conditions}. Skipping.")
        return None

    print(f"  Tumor group: {tumor_group}")
    print(f"  Normal group: {normal_group}")
    print(f"  Cancer cells: {adata_cancer.n_obs}")

    sc.tl.rank_genes_groups(
        adata_cancer,
        groupby="condition",
        groups=tumor_group,
        reference=normal_group[0],
        method="wilcoxon",
        use_raw=True,
        key_added="tumor_vs_normal",
    )

    result = sc.get.rank_genes_groups_df(adata_cancer, group=tumor_group[0],
                                          key="tumor_vs_normal")
    sig = result[result["pvals_adj"] < 0.05]
    print(f"\n  DEGs (FDR < 0.05): {len(sig):,}")

    return result


# ──────────────────────────────────────────────────────────────────────────────
# 4. Volcano plot
# ──────────────────────────────────────────────────────────────────────────────

def plot_volcano(
    de_result: pd.DataFrame,
    title: str = "Cancer Epithelial vs. Other Cell Types",
    filename: str = "deg_volcano.png",
    n_label: int = 15,
) -> None:
    """Create a publication-quality volcano plot."""
    print(f"\n[PLOT] Generating volcano plot: {filename} …")

    df = de_result.copy()
    df["neg_log10_padj"] = -np.log10(df["pvals_adj"].clip(lower=1e-300))

    # Classification
    df["significance"] = "Not significant"
    df.loc[
        (df["pvals_adj"] < 0.05) & (df["logfoldchanges"] > 1), "significance"
    ] = "Up in Cancer"
    df.loc[
        (df["pvals_adj"] < 0.05) & (df["logfoldchanges"] < -1), "significance"
    ] = "Down in Cancer"

    colors = {
        "Up in Cancer": "#e74c3c",
        "Down in Cancer": "#3498db",
        "Not significant": "#bdc3c7",
    }

    fig, ax = plt.subplots(figsize=(10, 8))

    for label, color in colors.items():
        subset = df[df["significance"] == label]
        ax.scatter(
            subset["logfoldchanges"],
            subset["neg_log10_padj"],
            c=color, alpha=0.5, s=8, label=label, edgecolors="none",
        )

    # Label top genes
    top_genes = pd.concat([
        df[df["significance"] == "Up in Cancer"].nlargest(n_label, "neg_log10_padj"),
        df[df["significance"] == "Down in Cancer"].nlargest(n_label // 2, "neg_log10_padj"),
    ])

    from matplotlib.offsetbox import AnnotationBbox

    for _, row in top_genes.iterrows():
        ax.annotate(
            row["names"],
            xy=(row["logfoldchanges"], row["neg_log10_padj"]),
            fontsize=7, fontweight="bold", alpha=0.8,
            xytext=(5, 5), textcoords="offset points",
        )

    ax.axhline(-np.log10(0.05), ls="--", color="grey", alpha=0.5, lw=0.8)
    ax.axvline(1, ls="--", color="grey", alpha=0.5, lw=0.8)
    ax.axvline(-1, ls="--", color="grey", alpha=0.5, lw=0.8)

    ax.set_xlabel("log₂ Fold Change", fontsize=12)
    ax.set_ylabel("-log₁₀ Adjusted p-value", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend(framealpha=0.9, fontsize=10)

    fig.tight_layout()
    fig.savefig(RESULTS_DIR / filename, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[PLOT] Saved → results/{filename}")


# ──────────────────────────────────────────────────────────────────────────────
# 5. Dot plot of top markers
# ──────────────────────────────────────────────────────────────────────────────

def plot_marker_dotplot(adata: ad.AnnData) -> None:
    """Dot plot of canonical markers across cell types."""
    print(f"\n[PLOT] Generating marker dot plot …")

    # Collect canonical markers that are present in the data
    # Import from the clustering module
    import importlib
    mod = importlib.import_module("scripts.02_clustering_and_annotation")
    MARKER_GENES = mod.MARKER_GENES

    var_names_in_data: dict[str, list[str]] = {}
    raw_names = set(adata.raw.var_names if adata.raw else adata.var_names)
    for ct, markers in MARKER_GENES.items():
        available = [m for m in markers if m in raw_names]
        if available:
            var_names_in_data[ct] = available[:3]  # top 3 per type

    if not var_names_in_data:
        print("  [WARN] No marker genes found in dataset for dot plot.")
        return

    fig = sc.pl.dotplot(
        adata, var_names=var_names_in_data, groupby="cell_type",
        use_raw=True, standard_scale="var", show=False, return_fig=True,
    )
    fig.savefig(RESULTS_DIR / "marker_dotplot.png", dpi=200, bbox_inches="tight")
    plt.close()
    print(f"[PLOT] Saved → results/marker_dotplot.png")


def plot_marker_dotplot_flat(adata: ad.AnnData) -> None:
    """
    Flat dot plot (not grouped) as fallback if the import-based approach fails.
    """
    # Canonical markers to show
    show_genes = [
        "EPCAM", "KRT19", "KRT8",       # Cancer epithelial
        "CD3D", "CD3E", "CD8A",          # T cells
        "CD79A", "MS4A1",                # B cells
        "CD68", "CD163",                 # Macrophages
        "COL1A1", "DCN",                 # Fibroblasts
        "PECAM1", "VWF",                 # Endothelial
        "NKG7", "GNLY",                  # NK cells
        "FCER1A", "CD1C",               # DCs
        "TPSAB1", "KIT",                # Mast cells
        "MZB1", "SDC1",                  # Plasmablasts
    ]

    raw_names = set(adata.raw.var_names if adata.raw else adata.var_names)
    available = [g for g in show_genes if g in raw_names]

    if not available:
        print("  [WARN] No marker genes found in dataset for dot plot.")
        return

    fig = sc.pl.dotplot(
        adata, var_names=available, groupby="cell_type",
        use_raw=True, standard_scale="var", show=False, return_fig=True,
    )
    fig.savefig(RESULTS_DIR / "marker_dotplot.png", dpi=200, bbox_inches="tight")
    plt.close()
    print(f"[PLOT] Saved → results/marker_dotplot.png")


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    if not INPUT_H5AD.exists():
        print(f"[ERROR] Input not found: {INPUT_H5AD}")
        print("        Run 02_clustering_and_annotation.py first.")
        return

    print(f"[INFO] Loading {INPUT_H5AD} …")
    adata = sc.read_h5ad(INPUT_H5AD)
    print(f"[INFO] Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    # ── Marker genes per cell type ────────────────────────────────────────
    markers_df = find_cell_type_markers(adata)

    # ── Cancer vs Other DE ────────────────────────────────────────────────
    cancer_de = cancer_vs_other_de(adata)
    if cancer_de is not None:
        plot_volcano(cancer_de)

    # ── Tumor vs Normal DE ────────────────────────────────────────────────
    tn_de = tumor_vs_normal_de(adata)
    if tn_de is not None:
        plot_volcano(
            tn_de,
            title="Tumor vs. Normal Adjacent (Cancer Epithelial)",
            filename="deg_volcano_tumor_vs_normal.png",
        )

    # ── Dot plot ──────────────────────────────────────────────────────────
    try:
        plot_marker_dotplot(adata)
    except Exception:
        plot_marker_dotplot_flat(adata)

    # ── Save ──────────────────────────────────────────────────────────────
    adata.write_h5ad(OUTPUT_H5AD)
    print(f"\n[INFO] Saved DE results → {OUTPUT_H5AD}")
    print("\n✅ Step 3 complete. Run 04_visualization_and_interpretation.py next.")


if __name__ == "__main__":
    main()
