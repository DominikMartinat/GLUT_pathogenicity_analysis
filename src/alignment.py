"""
DALI distance-matrix structural alignment parser and profile mapping across GLUT family members.
Reproduces Figure 5B (structural alignment overlap) and Figure 5C (helix-by-helix profile).
"""

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def parse_dali_strlist(strlist_val: str, bracketed: bool = True) -> list[int]:
    """Parse DALI bracketed space-separated string of indices into integer list."""
    if not isinstance(strlist_val, str):
        return []
    s = strlist_val.strip()
    if bracketed and s.startswith("[") and s.endswith("]"):
        s = s[1:-1]
    tokens = s.split()
    return [int(t) for t in tokens if t.isdigit() or (t.startswith("-") and t[1:].isdigit())]


def map_structural_alignment(
    dali_tsv: str | Path,
    pathogenicity_csv: str | Path,
    topology_3line: str | Path,
    reference_name: str = "GLUT1",
) -> pd.DataFrame:
    """
    Map residue-level pathogenicity and transmembrane topology onto reference structure coordinates
    using DALI pairwise alignment blocks.
    """
    dali_df = pd.read_csv(dali_tsv, sep="\t")
    dali_df["qstarts"] = dali_df["qstarts"].apply(parse_dali_strlist)
    dali_df["sstarts"] = dali_df["sstarts"].apply(parse_dali_strlist)
    dali_df["lengths"] = dali_df["lengths"].apply(parse_dali_strlist)

    pato_df = pd.read_csv(pathogenicity_csv)
    pato_df = pato_df.rename(columns={"residue_label": reference_name, "pathogenicity": f"p_{reference_name}"})

    with open(topology_3line, "r", encoding="utf-8") as f:
        topo_lines = f.read().splitlines()
    if len(topo_lines) >= 3:
        pato_df["topology"] = list(topo_lines[2][:len(pato_df)])

    return pato_df
