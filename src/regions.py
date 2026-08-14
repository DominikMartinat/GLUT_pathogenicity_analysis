"""
Region extraction and annotation tools for protein structures.
Supports:
- DeepTMHMM topology (intracellular 'I', extracellular 'O', membrane 'M')
- MOLEonline pore lining residues
- PrankWeb / P2Rank binding pocket residues
- All CA residues from PDB
"""

from pathlib import Path
import json
import pandas as pd


def extract_all_residues(pdb_file: str | Path, out_csv: str | Path) -> pd.DataFrame:
    """Extract all CA atom positions and residue names from a PDB structure."""
    residue_labels = []
    residue_names = []
    with open(pdb_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("ATOM"):
                parts = line.split()
                if len(parts) >= 6 and parts[2] == "CA":
                    residue_labels.append(int(parts[5]))
                    residue_names.append(parts[3])
    df = pd.DataFrame({"residue_label": residue_labels, "residue_name": residue_names})
    df.to_csv(out_csv, index=False)
    return df


def parse_deeptmhmm_topology(
    deeptmhmm_3line_file: str | Path,
    out_dir: str | Path,
    identifier: str,
    ext_prefix: str = "O",
    int_prefix: str = "I",
    mem_prefix: str = "M",
) -> tuple[Path, Path, Path]:
    """Parse DeepTMHMM 3line output and save separate CSVs for I, O, and M regions."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(deeptmhmm_3line_file, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()

    sequence = lines[1].strip()
    annotation = lines[2].strip()

    intracellular = {"residue_label": [], "residue_name": []}
    extracellular = {"residue_label": [], "residue_name": []}
    membrane = {"residue_label": [], "residue_name": []}

    for i, (res, reg) in enumerate(zip(sequence, annotation), start=1):
        if reg == "I":
            intracellular["residue_label"].append(i)
            intracellular["residue_name"].append(res)
        elif reg == "O":
            extracellular["residue_label"].append(i)
            extracellular["residue_name"].append(res)
        elif reg == "M":
            membrane["residue_label"].append(i)
            membrane["residue_name"].append(res)

    f_int = out_dir / f"{int_prefix}_{identifier}.csv"
    f_ext = out_dir / f"{ext_prefix}_{identifier}.csv"
    f_mem = out_dir / f"{mem_prefix}_{identifier}.csv"

    pd.DataFrame(intracellular).to_csv(f_int, index=False)
    pd.DataFrame(extracellular).to_csv(f_ext, index=False)
    pd.DataFrame(membrane).to_csv(f_mem, index=False)

    return f_int, f_ext, f_mem


def extract_binding_pocket_residues(prankweb_csv: str | Path, out_csv: str | Path) -> pd.DataFrame:
    """Extract primary (pocket #1) binding pocket residues from PrankWeb CSV output."""
    df = pd.read_csv(prankweb_csv)
    # Strip whitespace from column names
    df.columns = [c.strip() for c in df.columns]
    pocket1 = df[df["pocket"] == 1][["residue_label", "residue_name"]].copy()
    pocket1 = pocket1.sort_values(by="residue_label").drop_duplicates()
    pocket1.to_csv(out_csv, index=False)
    return pocket1


def extract_mole_lining_residues(moleonline_json: str | Path, out_csv: str | Path) -> pd.DataFrame:
    """Extract pore-lining residues from MOLEonline JSON output."""
    with open(moleonline_json, "r", encoding="utf-8") as f:
        data = json.load(f)

    lining_res_dict = {}
    paths = data.get("Channels", {}).get("Paths", [])
    if paths:
        residue_flow = paths[0].get("Layers", {}).get("ResidueFlow", [])
        for entry in residue_flow:
            tokens = entry.split()
            if len(tokens) >= 2:
                lining_res_dict[int(tokens[1])] = tokens[0]

    df = pd.DataFrame({
        "residue_label": list(lining_res_dict.keys()),
        "residue_name": list(lining_res_dict.values()),
    }).sort_values(by="residue_label")
    df.to_csv(out_csv, index=False)
    return df


def extract_lining_without_binding_pocket(
    lining_csv: str | Path,
    binding_pocket_csv: str | Path,
    out_csv: str | Path,
) -> pd.DataFrame:
    """Compute lining residues that do not overlap with the primary binding pocket."""
    lr_df = pd.read_csv(lining_csv)
    bp_df = pd.read_csv(binding_pocket_csv)

    bp_labels = set(bp_df["residue_label"].values)
    non_bp = lr_df[~lr_df["residue_label"].isin(bp_labels)].copy()
    non_bp = non_bp.sort_values(by="residue_label")
    non_bp.to_csv(out_csv, index=False)
    return non_bp
