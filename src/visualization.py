"""
Visualization utilities for GLUT pathogenicity profiles.
Generates publication-ready heatmaps and structural alignment profiles.
"""

from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def create_heatmap(
    input_file: str | Path,
    method_name: str,
    output_img: str | Path,
    inverse_color: bool = False,
    figsize: tuple[float, float] = (10.0, 6.0),
    dpi: int = 300,
) -> None:
    """
    Generate and save a publication-quality heatmap of regional pathogenicity scores.
    """
    df = pd.read_csv(input_file)
    if "identifier" in df.columns:
        df.set_index("identifier", inplace=True)
    elif "protein" in df.columns:
        df.set_index("protein", inplace=True)

    plt.figure(figsize=figsize)
    cmap = "coolwarm_r" if inverse_color else "coolwarm"

    sns.heatmap(
        df,
        cmap=cmap,
        vmin=0.0,
        vmax=1.0,
        linewidths=0.5,
        linecolor="gray",
        annot=True,
        fmt=".3f",
        cbar_kws={"label": f"{method_name} Score"},
    )

    plt.title(f"GLUT Family Pathogenicity Profile ({method_name})", fontsize=14, weight="bold")
    plt.xticks(rotation=45, ha="right", fontsize=11)
    plt.yticks(fontsize=11)
    plt.tight_layout()

    output_path = Path(output_img)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()
    print(f"Heatmap saved to {output_path}")
