"""
ClinVar benchmark script for GLUT pathogenicity analysis.
Loads curated ClinVar missense variants (Table 2), scores them against
AlphaMissense, PolyPhen-2, and SIFT, computes ROC-AUC curves (Figure 6B),
and plots the pathogenic-to-benign variant distribution across TM helices (Figure 6C).
"""

from pathlib import Path
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score


def run_clinvar_benchmark(data_csv: Path, output_dir: Path) -> dict[str, float]:
    """Execute ROC-AUC calculation and generate Figures 6B and 6C."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(data_csv)
    print(f"Loaded {len(df)} ClinVar variants across {df['protein'].nunique()} GLUT proteins.")
    print("Class breakdown:\n", df["clinical_significance"].value_counts())

    # Binary label: Pathogenic = 1, Benign = 0
    y_true = (df["clinical_significance"].str.lower().str.contains("pathogenic")).astype(int)

    # Predictions
    # Note: For SIFT, lower score = more deleterious, so inverted (1 - score) is used for standard ROC ranking
    scores = {
        "AlphaMissense": df["alphamissense_score"],
        "SIFT": 1.0 - df["sift_score"],
        "PolyPhen-2": df["polyphen2_score"],
    }

    aucs = {}
    curves = {}
    for name, y_score in scores.items():
        # Drop NaNs if any
        valid = y_score.notna()
        auc = roc_auc_score(y_true[valid], y_score[valid])
        fpr, tpr, _ = roc_curve(y_true[valid], y_score[valid])
        aucs[name] = auc
        curves[name] = (fpr, tpr)
        print(f"ROC-AUC ({name}): {auc:.4f}")

    # Plot Figure 6B: ROC Curves
    plt.figure(figsize=(7, 6))
    colors = {"AlphaMissense": "#d62728", "SIFT": "#2ca02c", "PolyPhen-2": "#1f77b4"}

    for name in ["AlphaMissense", "SIFT", "PolyPhen-2"]:
        fpr, tpr = curves[name]
        auc = aucs[name]
        plt.plot(fpr, tpr, label=f"{name} (AUC = {auc:.2f})", color=colors[name], lw=2.5)

    plt.plot([0, 1], [0, 1], "k--", lw=1.5, alpha=0.7, label="Chance (AUC = 0.50)")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate (1 - Specificity)", fontsize=12)
    plt.ylabel("True Positive Rate (Sensitivity)", fontsize=12)
    plt.title("Figure 6B: ROC Curves for Missense Pathogenicity Predictors", fontsize=13, weight="bold")
    plt.legend(loc="lower right", frameon=True, fontsize=11)
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()

    fig6b_path = output_dir / "figure6b_roc_curves.png"
    fig6b_svg = output_dir / "figure6b_roc_curves.svg"
    plt.savefig(fig6b_path, dpi=300, bbox_inches="tight")
    plt.savefig(fig6b_svg, bbox_inches="tight")
    plt.close()
    print(f"Saved Figure 6B to {fig6b_path} and {fig6b_svg}")

    # Plot Figure 6C: ClinVar variants per TM helix
    if "helix" in df.columns:
        helix_df = df[df["helix"].notna() & (df["helix"] != "loop") & (df["helix"] != "-")].copy()
        # Count pathogenic vs benign per helix
        helix_counts = (
            helix_df.groupby(["helix", "clinical_significance"])
            .size()
            .unstack(fill_value=0)
        )
        
        # Sort helices 1 to 12
        def helix_sort_key(h):
            try:
                digits = "".join([c for c in str(h) if c.isdigit()])
                return int(digits) if digits else 99
            except:
                return 99

        sorted_helices = sorted(helix_counts.index, key=helix_sort_key)
        helix_counts = helix_counts.reindex(sorted_helices).fillna(0)

        plt.figure(figsize=(10, 5))
        x = np.arange(len(sorted_helices))
        width = 0.38

        p_col = [c for c in helix_counts.columns if "pathogenic" in c.lower()]
        b_col = [c for c in helix_counts.columns if "benign" in c.lower()]

        p_vals = helix_counts[p_col[0]] if p_col else np.zeros(len(sorted_helices))
        b_vals = helix_counts[b_col[0]] if b_col else np.zeros(len(sorted_helices))

        plt.bar(x - width/2, p_vals, width, label="Pathogenic / Likely Pathogenic", color="#d62728", edgecolor="black")
        plt.bar(x + width/2, b_vals, width, label="Benign / Likely Benign", color="#1f77b4", edgecolor="black")

        plt.xlabel("Transmembrane Helix", fontsize=12)
        plt.ylabel("Number of ClinVar Variants", fontsize=12)
        plt.title("Figure 6C: ClinVar Variant Distribution Across Transmembrane Helices", fontsize=13, weight="bold")
        plt.xticks(x, [f"TM{h}" if not str(h).startswith("TM") else str(h) for h in sorted_helices], fontsize=11)
        plt.legend(loc="upper right", frameon=True, fontsize=11)
        plt.grid(axis="y", linestyle="--", alpha=0.5)
        plt.tight_layout()

        fig6c_path = output_dir / "figure6c_helix_variants.png"
        fig6c_svg = output_dir / "figure6c_helix_variants.svg"
        plt.savefig(fig6c_path, dpi=300, bbox_inches="tight")
        plt.savefig(fig6c_svg, bbox_inches="tight")
        plt.close()
        print(f"Saved Figure 6C to {fig6c_path} and {fig6c_svg}")

    return aucs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ClinVar benchmark and generate ROC curves")
    parser.add_argument("--data", type=Path, default=Path("data/clinvar/glut_clinvar_variants.csv"), help="Path to ClinVar variants CSV")
    parser.add_argument("--output-dir", type=Path, default=Path("results/figures"), help="Output directory for figures")
    args = parser.parse_args()

    run_clinvar_benchmark(args.data, args.output_dir)
