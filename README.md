# Mapping Pathogenic Patterns in Membrane Transporters from the GLUT Family

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.64898%2F2026.06.28.735151v1-red.svg)](https://www.biorxiv.org/content/10.64898/2026.06.28.735151v1)

This repository contains the complete computational workflow, structural models, regional mapping algorithms, and statistical analysis suite supporting the manuscript:

> **Mapping pathogenic patterns in membrane transporters from the GLUT transporter family**  
> Nina Kadášová, Dominik Martinát, Anna Špačková, Ivana Hutařová Vařeková, Karel Berka.  
> *bioRxiv* (2026). DOI: [10.64898/2026.06.28.735151v1](https://www.biorxiv.org/content/10.64898/2026.06.28.735151v1)

---

## Overview

The facilitative glucose transporter family (**GLUT**, encoded by *SLC2A1–SLC2A14*) plays an essential role in cellular energy metabolism. This project systematically maps missense mutation pathogenicity across:
1. **Topological Domains:** Transmembrane (TM) helices, intracellular loops, and extracellular segments (annotated via [DeepTMHMM](https://dtu.biolib.com/DeepTMHMM)).
2. **Translocation Pathways:** Pore-lining residues along the central cavity computed with [MOLEonline](https://mole.upol.cz).
3. **Substrate-Binding Sites:** Primary binding pockets predicted via [PrankWeb](https://prankweb.cz).
4. **Pathogenicity Predictors:** Benchmark comparisons between deep learning-based [AlphaMissense](https://github.com/google-deepmind/alphamissense), evolutionary-based [SIFT](https://sift.bii.a-star.edu.sg), and classifier-based [PolyPhen-2](http://genetics.bwh.harvard.edu/pph2) (via Rhapsody).
5. **Clinical Validation:** Empirical benchmark against curated human clinical missense variants from [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) with ROC-AUC analysis.

---

## Repository Structure

```text
GLUT_pathogenicity_analysis/
├── analyze_variant_impact.py      # Top-level CLI for regional statistics & Fig 4
├── requirements.txt               # Pip dependencies
├── environment.yml                # Conda environment definition
├── LICENSE                        # Creative Commons Attribution 4.0 (CC BY 4.0)
├── README.md                      # Documentation
├── data/
│   ├── statistics/                # Standardized regional pathogenicity TSVs (N=14)
│   │   ├── PyMissense.tsv         # AlphaMissense scores
│   │   ├── PolyPhen-2.tsv         # PolyPhen-2 scores
│   │   └── SIFT.tsv               # SIFT scores
│   ├── clinvar/                   # Curated ClinVar missense variants (Table 2)
│   │   └── glut_clinvar_variants.csv
│   └── data_backup/               # Raw inputs (PDBs, MOLE JSONs, PrankWeb CSVs, DALI TSVs)
├── src/                           # Modular Python package
│   ├── __init__.py
│   ├── analyze_variant_impact.py  # Statistical engine & paired t-test pipeline
│   ├── clinvar_benchmark.py       # ClinVar ROC-AUC & helix distribution module
│   ├── data_fetch.py              # FASTA/structure retrieval helpers
│   ├── regions.py                 # Topology & pore/pocket region parsing
│   ├── pathogenicity.py           # Residue-level pathogenicity assignment
│   ├── alignment.py               # DALI structural alignment parser
│   └── visualization.py           # Heatmap and plotting routines
├── notebooks/                     # Interactive, reproducible Jupyter notebooks
│   ├── 01_region_and_pathogenicity_mapping.ipynb   # Generates Figs 3 & 4
│   ├── 02_dali_structural_alignment.ipynb          # Generates Fig 5
│   └── 03_clinvar_clinical_benchmark.ipynb         # Generates Fig 6 & Table 2
└── results/
    └── figures/                   # Output figures (PNG & SVG) and summary statistics
        ├── figure4_region_comparison_bars.png
        ├── figure4_region_comparison_bars.svg
        ├── figure6b_roc_curves.png
        ├── figure6b_roc_curves.svg
        ├── figure6c_helix_variants.png
        ├── figure6c_helix_variants.svg
        ├── pvalue_heatmap_*.png
        └── region_summary_stats.csv
```

---

## Installation & Setup

### Option 1: Conda (Recommended)

```bash
git clone https://github.com/DominikMartinat/GLUT_pathogenicity_analysis.git
cd GLUT_pathogenicity_analysis
conda env create -f environment.yml
conda activate glut_pathogenicity
```

### Option 2: Virtual Environment (pip / uv)

```bash
python -m venv .venv
source .venv/bin/activate    # On Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

---

## Reproducing Paper Results

### 1. Statistical Analysis & Grouped Bar Chart (Figure 4)

Run the statistical pipeline on the 14 GLUT transporters ($N=14$):

```bash
python analyze_variant_impact.py --data-dir data/statistics --output-dir results/statistics
```

This generates:
* `results/figures/figure4_region_comparison_bars.png` (and `.svg`) — Grouped bar chart with paired $t$-test significance brackets ($^*p < 0.05, ^{**}p < 0.005, ^{***}p < 0.0005$).
* `results/figures/pvalue_heatmap_<method>.png` (and `.svg`) — Pairwise $p$-value heatmaps.
* `results/figures/region_summary_stats.csv` — Mean, SEM, and 95% Confidence Intervals per region.

### 2. Clinical Variant Benchmark & ROC Curves (Figure 6 & Table 2)

Run the ClinVar validation pipeline on the 140 curated clinical variants:

```bash
python src/clinvar_benchmark.py --data data/clinvar/glut_clinvar_variants.csv --output-dir results/clinvar_benchmark
```

This generates:
* `results/figures/figure6b_roc_curves.png` (and `.svg`) — ROC-AUC curves comparing AlphaMissense ($\text{AUC} = 0.88$), SIFT ($\text{AUC} = 0.79$), and PolyPhen-2 ($\text{AUC} = 0.78$).
* `results/figures/figure6c_helix_variants.png` (and `.svg`) — Pathogenic vs. benign variant counts across transmembrane helices TM1–12.

### 3. Interactive Notebooks

Launch Jupyter to explore the step-by-step pipeline:

```bash
jupyter notebook notebooks/
```

---

## Citation

If you use this repository or workflow in your research, please cite:

```bibtex
@article{kadasova2026glut,
  title={Mapping pathogenic patterns in membrane transporters from the GLUT transporter family},
  author={Kad{\'a}{\v{s}}ov{\'a}, Nina and Martin{\'a}t, Dominik and {\v{S}}pa{\v{c}}kov{\'a}, Anna and Huta{\v{r}}ov{\'a} Va{\v{r}}ekov{\'a}, Ivana and Berka, Karel},
  journal={bioRxiv},
  year={2026},
  doi={10.64898/2026.06.28.735151v1},
  url={https://www.biorxiv.org/content/10.64898/2026.06.28.735151v1}
}
```

---

## License

This project is licensed under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](LICENSE).