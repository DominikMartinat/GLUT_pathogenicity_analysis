"""
Pathogenicity extraction, mapping, and regional averaging for:
- AlphaMissense
- PolyPhen-2 (via Rhapsody)
- SIFT
"""

from pathlib import Path
import pandas as pd


def parse_alphamissense_pdb(am_pdb_file: str | Path, out_csv: str | Path) -> pd.DataFrame:
    """Extract residue-level AlphaMissense scores from an annotated PDB B-factor column."""
    residue_labels = []
    residue_names = []
    pathogenicities = []
    with open(am_pdb_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("ATOM"):
                parts = line.split()
                if len(parts) >= 11 and parts[2] == "CA":
                    residue_labels.append(int(parts[5]))
                    residue_names.append(parts[3])
                    pathogenicities.append(float(parts[10]))

    df = pd.DataFrame({
        "residue_label": residue_labels,
        "residue_name": residue_names,
        "pathogenicity": pathogenicities,
    }).drop_duplicates(subset=["residue_label"])
    df.to_csv(out_csv, index=False)
    return df


def parse_polyphen2_rhapsody(polyphen_file: str | Path, out_csv: str | Path) -> pd.DataFrame:
    """Parse PolyPhen-2 raw tabular output from Rhapsody and average across substitutions per site."""
    columns = [
        "ID", "residue_label", "residue_name", "AA_alt", "Training", "Score", "Prob", "Class",
        "pathogenicity", "PolyPhen_class", "EVmutation_score", "EVmutation_class"
    ]
    df = pd.read_csv(polyphen_file, sep=r"\s+", names=columns, comment="#", na_values="nan")
    df = df[["residue_label", "residue_name", "pathogenicity"]].dropna(subset=["pathogenicity"])

    result = df.groupby(["residue_label", "residue_name"]).agg({"pathogenicity": "mean"}).reset_index()
    result["pathogenicity"] = result["pathogenicity"].round(4)
    result.to_csv(out_csv, index=False)
    return result


def parse_sift_output(sift_file: str | Path, out_csv: str | Path) -> pd.DataFrame:
    """Parse SIFT matrix output and compute mean tolerated probability score per residue position."""
    names = "ACDEFGHIKLMNPQRSTVWY"
    active = False
    residue_labels = []
    residue_names = []
    pathogenicities = []

    with open(sift_file, "r", encoding="utf-8") as f:
        for line in f:
            line_str = line.strip()
            if not line_str:
                continue
            if line_str.startswith("pos"):
                active = True
                continue
            if active:
                parts = line_str.split()
                if len(parts) < 22:
                    continue
                pos_str = parts[0][:-1]
                native_aa = parts[0][-1]
                scores = parts[2:]
                non_native_sum = 0.0
                count = 0
                for i, aa in enumerate(names):
                    if aa != native_aa and i < len(scores):
                        non_native_sum += float(scores[i])
                        count += 1
                avg_score = non_native_sum / max(1, count)
                residue_labels.append(int(pos_str))
                residue_names.append(native_aa)
                pathogenicities.append(round(avg_score, 4))

    df = pd.DataFrame({
        "residue_label": residue_labels,
        "residue_name": residue_names,
        "pathogenicity": pathogenicities,
    })
    df.to_csv(out_csv, index=False)
    return df


def assign_pathogenicity_to_region(
    pathogenicity_csv: str | Path,
    region_csv: str | Path,
    out_csv: str | Path,
) -> pd.DataFrame:
    """Filter whole-protein pathogenicity scores to only include residues in a specific region."""
    patho_df = pd.read_csv(pathogenicity_csv)
    reg_df = pd.read_csv(region_csv)

    reg_labels = set(reg_df["residue_label"].values)
    filtered = patho_df[patho_df["residue_label"].isin(reg_labels)].copy()
    filtered.to_csv(out_csv, index=False)
    return filtered


def compute_average_pathogenicity(
    patho_dir: str | Path,
    region_prefix: str,
    id_list: list[str],
    out_csv: str | Path,
) -> pd.DataFrame:
    """Compute the mean pathogenicity for a given region across a list of protein identifiers."""
    patho_dir = Path(patho_dir)
    avg_scores = []
    identifiers = []
    for pid in id_list:
        file_path = patho_dir / f"{region_prefix}_{pid}.csv"
        if file_path.exists():
            df = pd.read_csv(file_path)
            avg_scores.append(df["pathogenicity"].mean())
        else:
            avg_scores.append(None)
        identifiers.append(pid)

    out_df = pd.DataFrame({"identifier": identifiers, "average_pathogenicity": avg_scores})
    out_df.to_csv(out_csv, index=False)
    return out_df
