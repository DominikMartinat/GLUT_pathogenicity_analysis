#!/usr/bin/env python
"""
analyze_variant_impact.py
==========================

Statistical analysis of missense-variant impact predictions (PolyPhen-2,
PyMissense / AlphaMissense, SIFT) across structural regions of the GLUT
membrane transporter family.

Reproduces:
- Figure 4: Grouped bar chart (Mean +/- 95% CI) with paired t-test significance brackets
- Supplementary Figure S1: P-value heatmaps per method
- Supplementary Table S1: Summary statistics (Mean, SEM, 95% CI, N)

Usage:
    python analyze_variant_impact.py --data-dir data/statistics --output-dir results/figures
"""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns
from matplotlib.patches import Patch

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DATASETS = [
    ("PolyPhen", "PolyPhen-2.tsv"),
    ("PyMissense", "PyMissense.tsv"),
    ("SIFT", "SIFT.tsv"),
]

REGION_ORDER = [
    "binding place",
    "lining residues",
    "lining residues without binding place",
    "average for protein",
    "transmembrane region",
    "intracellular domain",
    "extracellular domain",
]

BINDING_SITE_GROUP = [
    "binding place",
    "lining residues",
    "lining residues without binding place",
    "average for protein",
]
TOPOLOGY_GROUP = [
    "transmembrane region",
    "average for protein",
    "intracellular domain",
    "extracellular domain",
]

BAR_WIDTH = 0.15
GROUP_GAP = 0.5
SIGNIFICANCE_BINS = [0, 0.0005, 0.005, 0.05, 1]
SIGNIFICANCE_LABELS = ["***", "**", "*", "ns"]


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_dataset(filepath: Path, dataset_name: str) -> pd.DataFrame:
    """Load one tab-separated region-score file and tag it with its method name."""
    df = pd.read_csv(filepath, sep="\t")
    df["Dataset"] = dataset_name
    return df


def load_all_datasets(data_dir: Path) -> dict[str, pd.DataFrame]:
    """Load every dataset listed in DATASETS from data_dir."""
    datasets = {}
    for name, filename in DATASETS:
        path = data_dir / filename
        if not path.exists():
            raise FileNotFoundError(
                f"Expected input file not found: {path}\n"
                f"See the module docstring / README.md for the required format."
            )
        datasets[name] = load_dataset(path, name)
    return datasets


def to_long_format(datasets: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Concatenate all datasets and reshape to long format (one row per protein x region)."""
    combined = pd.concat(datasets.values(), ignore_index=True)
    long_df = pd.melt(
        combined,
        id_vars=["protein", "Dataset"],
        var_name="Region",
        value_name="Value",
    )
    return long_df


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def compute_group_stats(long_df: pd.DataFrame, confidence: float = 0.95) -> pd.DataFrame:
    """Compute the mean and confidence interval of Value for every Dataset x Region group."""
    results = []
    for (dataset, region), group in long_df.groupby(["Dataset", "Region"]):
        values = group["Value"].dropna()
        mean = values.mean()
        sem = stats.sem(values)
        ci = sem * stats.t.ppf((1 + confidence) / 2.0, len(values) - 1)
        results.append({
            "Dataset": dataset,
            "Region": region,
            "Mean": mean,
            "SEM": sem,
            "CI": ci,
            "N": len(values),
        })

    summary_df = pd.DataFrame(results)
    summary_df["Region"] = pd.Categorical(summary_df["Region"], categories=REGION_ORDER, ordered=True)
    return summary_df.sort_values(["Dataset", "Region"]).reset_index(drop=True)


def compute_pvalue_matrix(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """Paired t-test p-values between every pair of region columns (same protein, two-tailed)."""
    pval_matrix = pd.DataFrame(np.nan, index=columns, columns=columns)
    for col1, col2 in combinations(columns, 2):
        if col1 in df.columns and col2 in df.columns:
            _, p = stats.ttest_rel(df[col1], df[col2])
            pval_matrix.loc[col1, col2] = p
            pval_matrix.loc[col2, col1] = p
    return pval_matrix


def pvalue_to_stars(p: float) -> str:
    """Convert a p-value to a conventional significance annotation."""
    if p is None or pd.isna(p):
        return ""
    if p < 0.0005:
        return "***"
    if p < 0.005:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def pvalue_matrix_to_dict(pval_matrix: pd.DataFrame) -> dict[tuple[str, str], float]:
    """Flatten a symmetric p-value matrix into a {sorted_pair: p_value} lookup."""
    pval_dict = {}
    for col1, col2 in combinations(pval_matrix.columns, 2):
        pval_dict[tuple(sorted((col1, col2)))] = pval_matrix.loc[col1, col2]
    return pval_dict


# ---------------------------------------------------------------------------
# Plot 1: Grouped bar chart with significance brackets (Figure 4)
# ---------------------------------------------------------------------------

def _calc_positions(start: float, regions: list[str]) -> list[float]:
    return [start + i * BAR_WIDTH for i in range(len(regions))]


def _build_group_positions(dataset_names: list[str]) -> list[list[float]]:
    """X positions for each dataset's two region sub-groups, with a small gap between them."""
    group_positions = []
    start = 0.0
    for _ in dataset_names:
        pos_binding = _calc_positions(start, BINDING_SITE_GROUP)
        start = pos_binding[-1] + (BAR_WIDTH * 2)
        pos_topology = _calc_positions(start, TOPOLOGY_GROUP)
        start = pos_topology[-1] + GROUP_GAP
        group_positions.append(pos_binding + pos_topology)
    return group_positions


def _draw_bracket_with_pvalue(ax, x1, x2, y, text, height=0.03, color="black", lw=0.8):
    """Draw a bracket between x1 and x2 at height y, with a significance label above it."""
    ax.plot([x1, x1, x2, x2], [y + 0.005, y + 0.025, y + 0.025, y + 0.005], color=color, lw=lw)
    ax.text((x1 + x2) / 2, y + height, text, ha="center", va="bottom", fontsize=14, color=color, weight="bold")


def plot_grouped_bars_with_significance(
    summary_df: pd.DataFrame,
    pvalue_dicts: dict[str, dict[tuple[str, str], float]],
    output_path: Path | None = None,
    show: bool = False,
) -> None:
    """Grouped bar chart (mean +/- 95% CI) per dataset, annotated with pairwise significance brackets."""
    all_regions = BINDING_SITE_GROUP + TOPOLOGY_GROUP
    region_order_unique = list(dict.fromkeys(all_regions))

    tab10 = plt.get_cmap("tab10").colors
    region_color_map = {}
    for i, region in enumerate(region_order_unique):
        if region == "average for protein":
            region_color_map[region] = ("#d62728", "#8c1b1b")
        else:
            region_color_map[region] = (tab10[i % len(tab10)], "black")

    dataset_names = list(summary_df["Dataset"].unique())
    group_positions = _build_group_positions(dataset_names)

    fig, ax = plt.subplots(figsize=(15, 7.5))

    for d_idx, dataset in enumerate(dataset_names):
        positions = group_positions[d_idx]

        for r_idx, x_pos in enumerate(positions):
            region = all_regions[r_idx]
            row = summary_df[(summary_df["Dataset"] == dataset) & (summary_df["Region"] == region)]
            if row.empty:
                continue
            mean, ci = row["Mean"].values[0], row["CI"].values[0]
            bar_color, edge_color = region_color_map[region]
            is_average = region == "average for protein"
            ax.bar(
                x_pos, mean, width=BAR_WIDTH, yerr=ci, capsize=4,
                label=region if d_idx == 0 and r_idx < len(BINDING_SITE_GROUP) else "",
                color=bar_color, edgecolor=edge_color,
                alpha=1.0 if is_average else 0.85,
                linewidth=1.8 if is_average else 1.0,
            )

        # Significance brackets
        pval_dict = pvalue_dicts.get(dataset, {})
        for i in range(len(positions) - 1):
            region1, region2 = all_regions[i], all_regions[i + 1]
            if region1 == BINDING_SITE_GROUP[-1] and region2 == TOPOLOGY_GROUP[0]:
                continue
            row1 = summary_df[(summary_df["Dataset"] == dataset) & (summary_df["Region"] == region1)]
            row2 = summary_df[(summary_df["Dataset"] == dataset) & (summary_df["Region"] == region2)]
            if row1.empty or row2.empty:
                continue
            y = max(row1["Mean"].values[0] + row1["CI"].values[0], row2["Mean"].values[0] + row2["CI"].values[0]) + 0.04
            p_val = pval_dict.get(tuple(sorted((region1, region2))))
            stars = pvalue_to_stars(p_val)
            if stars != "ns":
                _draw_bracket_with_pvalue(ax, positions[i], positions[i + 1], y, stars)

    xtick_positions = [(pos[0] + pos[-1]) / 2 for pos in group_positions]
    ax.set_xticks(xtick_positions)
    ax.set_xticklabels(dataset_names, fontsize=15, weight="bold")
    ax.set_ylabel("Mean Pathogenicity Score ± 95% CI", fontsize=14)
    ax.set_title("Figure 4: Regional Pathogenicity Comparison Across GLUT Transporters (N=14)", fontsize=15, weight="bold")
    ax.tick_params(axis="both", labelsize=12)

    legend_patches = [
        Patch(facecolor=region_color_map[r][0], edgecolor=region_color_map[r][1], label=r.title())
        for r in region_order_unique
    ]
    ax.legend(handles=legend_patches, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True, fontsize=11)
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    fig.tight_layout()

    if output_path:
        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_p, dpi=300, bbox_inches="tight")
        fig.savefig(out_p.with_suffix(".svg"), bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


# ---------------------------------------------------------------------------
# Plot 2: P-value Heatmap
# ---------------------------------------------------------------------------

def plot_pvalue_heatmap(pval_matrix: pd.DataFrame, title: str, output_path: Path | None = None, show: bool = False) -> None:
    """Heatmap of pairwise p-values, colored by significance band and annotated with stars."""
    colors = ["#1a9850", "#fee08b", "#fc8d59", "#d73027"]
    cmap = mcolors.ListedColormap(colors)
    norm = mcolors.BoundaryNorm(SIGNIFICANCE_BINS, cmap.N)

    annotations = pval_matrix.copy().astype(object)
    for i in annotations.index:
        for j in annotations.columns:
            p = pval_matrix.loc[i, j]
            annotations.loc[i, j] = f"{p:.1e}\n{pvalue_to_stars(p)}" if pd.notna(p) else ""

    fig, ax = plt.subplots(figsize=(9, 7.5))
    sns.heatmap(
        pval_matrix, cmap=cmap, norm=norm, annot=annotations, fmt="",
        linewidths=0.5, cbar_kws={"label": "P-value"}, square=True,
        mask=pval_matrix.isnull(), ax=ax,
    )
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks(SIGNIFICANCE_BINS)
    colorbar.set_ticklabels(["0", "0.0005", "0.005", "0.05", "1"])

    ax.set_title(title, fontsize=13, weight="bold")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=10)
    fig.tight_layout()

    if output_path:
        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_p, dpi=300, bbox_inches="tight")
        fig.savefig(out_p.with_suffix(".svg"), bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(data_dir: Path, output_dir: Path, show: bool) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    datasets = load_all_datasets(data_dir)
    long_df = to_long_format(datasets)
    summary_df = compute_group_stats(long_df)

    pvalue_dicts = {}
    for dataset_name, df in datasets.items():
        pval_matrix = compute_pvalue_matrix(df, REGION_ORDER)
        pvalue_dicts[dataset_name] = pvalue_matrix_to_dict(pval_matrix)
        plot_pvalue_heatmap(
            pval_matrix,
            title=f"Supplementary Fig S1: Pairwise Paired t-test P-values ({dataset_name})",
            output_path=output_dir / f"pvalue_heatmap_{dataset_name}.png",
            show=show,
        )

    plot_grouped_bars_with_significance(
        summary_df,
        pvalue_dicts,
        output_path=output_dir / "figure4_region_comparison_bars.png",
        show=show,
    )

    summary_df.to_csv(output_dir / "region_summary_stats.csv", index=False)
    print(f"Done. Figure 4, heatmaps, and summary stats written to {output_dir}/")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-dir", type=Path, default=Path("data/statistics"), help="Directory containing input .tsv files")
    parser.add_argument("--output-dir", type=Path, default=Path("results/figures"), help="Output directory for figures")
    parser.add_argument("--show", action="store_true", help="Display figures interactively")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args.data_dir, args.output_dir, args.show)
