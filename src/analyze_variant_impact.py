#!/usr/bin/env python
"""
analyze_variant_impact.py
==========================

Statistical analysis of missense-variant impact predictions (PolyPhen-2,
AlphaMissense / PyMissense, SIFT) across structural regions of the GLUT
membrane transporter family (SLC2A1–SLC2A14).

Features:
1. Robust data loading with auto-detection of delimiters and decimal separators (. or ,).
2. Canonical header aliasing (e.g. 'all' -> 'average for protein', 'binding pocket' -> 'binding place').
3. Paired t-tests between neighboring structural regions along the functional/topological gradient:
   binding place -> lining residues -> transmembrane region -> average for protein -> intracellular -> extracellular.
4. Publication-grade grouped bar chart (Mean +/- 95% CI) with pairwise significance brackets (Figure 4).
5. Full pairwise p-value heatmaps per method (Supplementary Figure S1).
6. Comprehensive export of summary statistics and p-values to CSV and Excel.

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
    ("PolyPhen", "PolyPhen-2.tsv", "PolyPhen2_NEW.txt"),
    ("AlphaMissense", "PyMissense.tsv", "AlphaMissense_NEW.txt"),
    ("SIFT", "SIFT.tsv", "SIFT_NEW.txt"),
]

REGIONS_ORDERED = [
    "binding place",
    "lining residues",
    "transmembrane region",
    "average for protein",
    "intracellular domain",
    "extracellular domain",
]

BAR_WIDTH = 0.14
GAP_BETWEEN_DATASETS = 0.45
SHOW_NS = True

SIGNIFICANCE_BINS = [0, 0.0005, 0.005, 0.05, 1]
SIGNIFICANCE_LABELS = ["***", "**", "*", "ns"]

COLUMN_ALIASES = {
    "protein": "protein",
    "identifier": "protein",
    "average for protein": "average for protein",
    "all": "average for protein",
    "binding place": "binding place",
    "binding pocket": "binding place",
    "lining residues": "lining residues",
    "lining residues without binding place": "lining residues without binding place",
    "transmembrane region": "transmembrane region",
    "transmembrane": "transmembrane region",
    "intracellular domain": "intracellular domain",
    "intracellular": "intracellular domain",
    "extracellular domain": "extracellular domain",
    "extracellular": "extracellular domain",
}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _read_table(filepath: Path, decimal: str) -> pd.DataFrame:
    """Read a table file, falling back to Python engine if needed."""
    try:
        df = pd.read_csv(filepath, sep="\t", decimal=decimal)
        if df.shape[1] == 1:
            df = pd.read_csv(filepath, sep=None, engine="python", decimal=decimal)
    except Exception:
        df = pd.read_csv(filepath, sep=None, engine="python", decimal=decimal)
    df.columns = df.columns.astype(str).str.strip()
    return df


def _rename_to_canonical_columns(df: pd.DataFrame, filepath: Path) -> pd.DataFrame:
    """Rename recognized header variants to canonical column names."""
    rename_map = {}
    for col in df.columns:
        key = str(col).strip().lower()
        if key in COLUMN_ALIASES:
            rename_map[col] = COLUMN_ALIASES[key]
    df = df.rename(columns=rename_map)
    if "protein" not in df.columns:
        raise ValueError(
            f"{filepath}: could not find a protein identifier column "
            f"(expected 'protein' or 'identifier'). Found: {list(df.columns)}"
        )
    return df


def load_one_file(filepath: Path, dataset_name: str) -> pd.DataFrame:
    """Load one region-score file, handling either . or , decimal format."""
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    df = _read_table(filepath, decimal=".")
    df = _rename_to_canonical_columns(df, filepath)

    region_columns = [c for c in df.columns if c != "protein"]
    if any(not pd.api.types.is_numeric_dtype(df[c]) for c in region_columns):
        df = _read_table(filepath, decimal=",")
        df = _rename_to_canonical_columns(df, filepath)

    # Cast all region columns to float
    for c in region_columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df["Dataset"] = dataset_name
    return df


def load_all_datasets(data_dir: Path) -> dict[str, pd.DataFrame]:
    """Load every dataset listed in DATASETS from data_dir."""
    datasets = {}
    for name, clean_file, raw_file in DATASETS:
        path = data_dir / clean_file
        if not path.exists():
            path = data_dir / raw_file
        if not path.exists():
            raise FileNotFoundError(f"Neither {clean_file} nor {raw_file} found in {data_dir}")
        datasets[name] = load_one_file(path, name)
    return datasets


def to_long_format(datasets: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Concatenate all datasets and reshape to long format."""
    combined = pd.concat(datasets.values(), ignore_index=True)
    long_df = combined.melt(
        id_vars=["protein", "Dataset"],
        var_name="Region",
        value_name="Value",
    )
    long_df["Region"] = long_df["Region"].astype(str).str.strip()
    long_df["Value"] = pd.to_numeric(long_df["Value"], errors="coerce")
    long_df = long_df.dropna(subset=["Value"])
    long_df = long_df[long_df["Region"].isin(REGIONS_ORDERED)].copy()
    return long_df


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def compute_group_stats(long_df: pd.DataFrame, confidence: float = 0.95) -> pd.DataFrame:
    """Compute the mean, SEM, and 95% confidence interval for each Dataset x Region group."""
    results = []
    for (dataset, region), group in long_df.groupby(["Dataset", "Region"]):
        values = group["Value"].dropna()
        n = len(values)
        mean = values.mean()
        sem = stats.sem(values) if n > 1 else 0.0
        ci = sem * stats.t.ppf((1 + confidence) / 2.0, n - 1) if n > 1 else 0.0
        results.append({
            "Dataset": dataset,
            "Region": region,
            "Mean": mean,
            "SEM": sem,
            "CI": ci,
            "N": n,
        })

    summary_df = pd.DataFrame(results)
    summary_df["Region"] = pd.Categorical(summary_df["Region"], categories=REGIONS_ORDERED, ordered=True)
    return summary_df.sort_values(["Dataset", "Region"]).reset_index(drop=True)


def p_to_stars(p: float) -> str:
    """Convert a p-value to standard significance annotation."""
    if p is None or pd.isna(p):
        return ""
    if p < 0.0005:
        return "***"
    if p < 0.005:
        return "**"
    if p < 0.05:
        return "*"
    return "ns" if SHOW_NS else ""


def compute_neighbor_pvalues(long_df: pd.DataFrame, datasets: list[str], regions_ordered: list[str]) -> pd.DataFrame:
    """Compute paired t-tests between adjacent/neighboring structural regions."""
    results = []
    for dataset in datasets:
        df_d = long_df[long_df["Dataset"] == dataset]
        for i in range(len(regions_ordered) - 1):
            r1 = regions_ordered[i]
            r2 = regions_ordered[i + 1]

            v1 = df_d[df_d["Region"] == r1].set_index("protein")["Value"]
            v2 = df_d[df_d["Region"] == r2].set_index("protein")["Value"]
            paired = pd.concat([v1, v2], axis=1, keys=["v1", "v2"]).dropna()

            if len(paired) > 1:
                stat, p_val = stats.ttest_rel(paired["v1"], paired["v2"])
            else:
                p_val = np.nan

            results.append({
                "Dataset": dataset,
                "Region_1": r1,
                "Region_2": r2,
                "N": len(paired),
                "p_value": p_val,
                "significance": p_to_stars(p_val),
            })
    return pd.DataFrame(results)


def compute_full_pvalue_matrix(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """Compute pairwise paired t-test p-values for all region pairs."""
    pval_matrix = pd.DataFrame(np.nan, index=columns, columns=columns)
    for col1, col2 in combinations(columns, 2):
        if col1 in df.columns and col2 in df.columns:
            paired = df[[col1, col2]].dropna()
            if len(paired) > 1:
                _, p = stats.ttest_rel(paired[col1], paired[col2])
                pval_matrix.loc[col1, col2] = p
                pval_matrix.loc[col2, col1] = p
    return pval_matrix


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def draw_bracket(ax, x1: float, x2: float, y: float, h: float, text: str, fontsize: int = 12):
    """Draw a statistical significance bracket between x1 and x2."""
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="black", linewidth=1.2)
    ax.text((x1 + x2) / 2, y + h + 0.005, text, ha="center", va="bottom", fontsize=fontsize, fontweight="bold")


def plot_grouped_bars_with_neighbor_significance(
    summary_df: pd.DataFrame,
    neighbor_pvals_df: pd.DataFrame,
    output_path: Path | None = None,
    show: bool = False,
) -> None:
    """Generate Figure 4: Grouped bar chart (Mean +/- 95% CI) with neighboring significance brackets."""
    tab10 = plt.get_cmap("tab10").colors
    region_color_map = {}
    for i, region in enumerate(REGIONS_ORDERED):
        if region == "average for protein":
            region_color_map[region] = ("#d62728", "#8c1b1b")
        else:
            region_color_map[region] = (tab10[i % len(tab10)], "black")

    datasets = list(summary_df["Dataset"].unique())

    dataset_positions = []
    start = 0.0
    for dataset in datasets:
        positions_for_dataset = []
        for region in REGIONS_ORDERED:
            positions_for_dataset.append((start, region))
            start += BAR_WIDTH
        dataset_positions.append(positions_for_dataset)
        start += GAP_BETWEEN_DATASETS

    fig, ax = plt.subplots(figsize=(15, 7.5))

    bar_tops = {}
    bar_x = {}
    max_y = 0.0

    for d_idx, (dataset, positions) in enumerate(zip(datasets, dataset_positions)):
        for x_pos, region in positions:
            row = summary_df[(summary_df["Dataset"] == dataset) & (summary_df["Region"] == region)]
            if row.empty:
                continue

            mean = row["Mean"].iloc[0]
            ci = row["CI"].iloc[0]

            bar_tops[(dataset, region)] = mean + ci
            bar_x[(dataset, region)] = x_pos
            max_y = max(max_y, mean + ci)

            bar_color, edge_color = region_color_map[region]
            is_average = region == "average for protein"

            ax.bar(
                x_pos,
                mean,
                width=BAR_WIDTH,
                yerr=ci,
                capsize=4,
                color=bar_color,
                edgecolor=edge_color,
                alpha=1.0 if is_average else 0.85,
                linewidth=1.8 if is_average else 1.0,
            )

    # Brackets
    bracket_gap = max(0.02, max_y * 0.035)
    bracket_height = max(0.015, max_y * 0.02)

    for _, row in neighbor_pvals_df.iterrows():
        dataset = row["Dataset"]
        r1 = row["Region_1"]
        r2 = row["Region_2"]
        sig = row["significance"]

        if not sig:
            continue

        x1 = bar_x.get((dataset, r1))
        x2 = bar_x.get((dataset, r2))
        y1 = bar_tops.get((dataset, r1))
        y2 = bar_tops.get((dataset, r2))

        if x1 is None or x2 is None:
            continue

        y = max(y1, y2) + bracket_gap
        draw_bracket(ax, x1, x2, y, bracket_height, sig, fontsize=12)

    # X-axis ticks
    xtick_positions = [(pos[0][0] + pos[-1][0]) / 2 for pos in dataset_positions]
    ax.set_xticks(xtick_positions)
    ax.set_xticklabels(datasets, fontsize=15, weight="bold")
    ax.set_ylabel("Mean Pathogenicity Score ± 95% CI", fontsize=14)
    ax.set_title("Figure 4: Regional Pathogenicity Comparison Across GLUT Transporters (N=14)", fontsize=15, weight="bold")
    ax.tick_params(axis="both", labelsize=12)

    legend_patches = [
        Patch(facecolor=region_color_map[r][0], edgecolor=region_color_map[r][1], label=r.title())
        for r in REGIONS_ORDERED
    ]
    ax.legend(handles=legend_patches, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True, fontsize=11)
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    ax.set_ylim(0, max_y * 1.25)
    fig.tight_layout()

    if output_path:
        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_p, dpi=300, bbox_inches="tight")
        fig.savefig(out_p.with_suffix(".svg"), bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def plot_pvalue_heatmap(pval_matrix: pd.DataFrame, title: str, output_path: Path | None = None, show: bool = False) -> None:
    """Supplementary Figure S1: Heatmap of pairwise p-values per predictor."""
    colors = ["#1a9850", "#fee08b", "#fc8d59", "#d73027"]
    cmap = mcolors.ListedColormap(colors)
    norm = mcolors.BoundaryNorm(SIGNIFICANCE_BINS, cmap.N)

    annotations = pval_matrix.copy().astype(object)
    for i in annotations.index:
        for j in annotations.columns:
            p = pval_matrix.loc[i, j]
            annotations.loc[i, j] = f"{p:.1e}\n{p_to_stars(p)}" if pd.notna(p) else ""

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
# Main Runner
# ---------------------------------------------------------------------------

def run(data_dir: Path, output_dir: Path, show: bool = False) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    datasets = load_all_datasets(data_dir)
    long_df = to_long_format(datasets)
    summary_df = compute_group_stats(long_df)
    dataset_names = list(datasets.keys())

    # 1. Neighboring comparisons & Figure 4
    neighbor_pvals = compute_neighbor_pvalues(long_df, dataset_names, REGIONS_ORDERED)
    plot_grouped_bars_with_neighbor_significance(
        summary_df,
        neighbor_pvals,
        output_path=output_dir / "figure4_region_comparison_bars.png",
        show=show,
    )

    # 2. Pairwise p-value heatmaps (Supplementary Figure S1)
    for name, df in datasets.items():
        pval_mat = compute_full_pvalue_matrix(df, REGIONS_ORDERED)
        plot_pvalue_heatmap(
            pval_mat,
            title=f"Supplementary Fig S1: Pairwise Paired t-test P-values ({name})",
            output_path=output_dir / f"pvalue_heatmap_{name}.png",
            show=show,
        )

    # 3. Export tables in CSV and Excel
    summary_df.to_csv(output_dir / "region_summary_stats.csv", index=False)
    summary_df.to_excel(output_dir / "summary_statistics_with_CI.xlsx", index=False)
    neighbor_pvals.to_csv(output_dir / "neighbor_region_pvalues.csv", index=False)
    neighbor_pvals.to_excel(output_dir / "neighbor_region_pvalues.xlsx", index=False)

    print(f"Analysis complete. Figures and tables generated in: {output_dir}/")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-dir", type=Path, default=Path("data/statistics"), help="Path to data directory")
    parser.add_argument("--output-dir", type=Path, default=Path("results/figures"), help="Path to output directory")
    parser.add_argument("--show", action="store_true", help="Display figures interactively")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args.data_dir, args.output_dir, args.show)
