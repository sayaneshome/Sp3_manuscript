# Sp3 RNA-seq analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Reproducible code for the Sp3 bulk RNA-seq study — **differential expression**,
**transcription-factor motif enrichment**, **GO/GSEA**, **figure generation**, and a
**machine-learning** analysis of TF-target gene expression. Comparisons are
Sp3-null vs. control in two cell types: **MEF** (mouse embryonic fibroblasts) and
**LMC**. Three contrasts are referred to throughout as **I, II, III**; "up"/"down"
denote up-/down-regulated gene sets in the Sp3 mutant.

This repository contains the **canonical final code only**, the **processed input
data**, and a **pinned Python + R environment**, so the analyses and figures can be
reproduced end-to-end.

---

## Quick start

```bash
# 1. Clone
git clone https://github.com/sayaneshome/sp3-rnaseq-analysis.git
cd sp3-rnaseq-analysis

# 2. Build the combined Python + R environment (one step)
conda env create -f environment/environment.yml   # or: mamba env create -f ...
conda activate sp3-rnaseq

# 3. Run — ALWAYS use data/ as the working directory
cd data
jupyter lab        # open ../bioinformatics/... or ../machine_learning/... notebooks
# or for an R script:
Rscript ../bioinformatics/plots_and_figures/mef_plots.R
```

> **Important:** all paths in the code are relative to the `data/` directory, so
> launch Jupyter (and run R scripts) **from inside `data/`**. See `data/README.md`.

Prefer pip + a system R? Use `environment/requirements.txt` (Python) and
`environment/install_R_packages.R` (R/Bioconductor) instead of the conda file.

---

## Repository layout

```
sp3-rnaseq-analysis/
├── README.md
├── LICENSE                       # MIT
├── CITATION.cff
├── CODE_AVAILABILITY.md          # manuscript-ready statement
├── environment/
│   ├── environment.yml           # one-step conda env (Python + R + Bioconductor)
│   ├── requirements.txt          # Python deps (pinned)
│   ├── install_R_packages.R      # CRAN + Bioconductor installer
│   └── R_sessionInfo_notes.md
├── data/                         # processed inputs (run code from here)
│   ├── README.md
│   └── _MANIFEST.txt
├── bioinformatics/
│   ├── differential_expression/  # DESeq2 DE + manuscript figure notebooks
│   ├── motif_analysis/           # HOMER motif enrichment + TF binding sites
│   └── plots_and_figures/        # GO/GSEA, scatterplots, Venn, heatmaps
├── machine_learning/             # random-forest model (manuscript Suppl. Fig 2c)
└── supplementary_tables/         # DE genes + GSEA pathway tables (per comparison)
```

**Supplementary tables** ([`supplementary_tables/`](supplementary_tables/)):
Table S1 — differentially expressed genes (gene, log2FC, p-value, adjusted p-value;
DESeq2, padj ≤ 0.05) and Table S2 — clusterProfiler GSEA pathway enrichment, each
with one sheet per comparison (LMC & MEF, I/II/III). Provided as CSVs and Excel
workbooks.

External command-line tools used upstream (their outputs are bundled in `data/`):
**HOMER** (motif enrichment) and **MEME Suite / FIMO** (motif scanning).

---

## Analysis pipeline & run order

Stages are sequential; each consumes the previous stage's outputs.

```
[1] Differential expression (DESeq2)  ──►  normalized counts, DE gene lists (I/II/III × up/down)
[2] Motif analysis (HOMER / FIMO)      ──►  per-contrast motif tables, TF-family summaries
[3] Figures & enrichment (GO/GSEA, scatterplots, Venn, heatmaps)
[4] Machine learning (random forest on TF-family features)  ──►  Supplementary Figure 2c
```

### Stage 1 — Differential expression (`bioinformatics/differential_expression/`)
| Notebook | Role |
|---|---|
| `2025_MEF_final_manuscript_svg_and_full_analyses.ipynb` | **Master MEF analysis** + manuscript figure SVGs |
| `2025_LMC_final_manuscript_svg_and_full_analyses.ipynb` | **Master LMC analysis** + manuscript figure SVGs |
| `DE_analysis_2023_MEF_2020counts_final_ongoing_with_svg_ggsave.ipynb` | Standalone MEF DESeq2 (PCA, volcano, heatmaps) |
| `DE_analysis_2023_LMC_2020counts_final_ongoing.ipynb` | Standalone LMC DESeq2 |

The `2025_*_final_manuscript_*` notebooks are the canonical entry points (DE + figures).
DESeq2 runs in R inside these notebooks via `rpy2` `%%R` cell magic.

### Stage 2 — Motif analysis (`bioinformatics/motif_analysis/`)
`fasta_sequence_derive.ipynb` (promoter sequences) → `Motif_scanning.ipynb`,
`MEF_homer.ipynb` / `LMC_homer.ipynb` (summarize HOMER enrichment) →
`homer_plots_I_II_analysis.R` / `homer_plots_III_analysis.R` (TF-family plots),
`TF_binding_sites_code.R` (binding-site annotation).

### Stage 3 — Figures & enrichment (`bioinformatics/plots_and_figures/`)
`clusterprofiler_GO_sorting_figure2_and_3.R` (GO/GSEA, Figs 2 & 3),
`code_LMC_MEF_normalised_difference_scatterplots*.R` (scatterplots),
`GC_TATA_venndiagrams_MEF_LMC.ipynb` (GC/TATA Venn), `heatmap_version_merging.ipynb`,
`mef_plots.R`.

### Stage 4 — Machine learning (`machine_learning/`)
**`ML_analysis_manuscript_Feb17_2-26.ipynb`** — the manuscript machine-learning
analysis, producing **Supplementary Figure 2c**. A random-forest regressor predicts
TF-target gene log2 fold-change from transcription-factor–family features, evaluated
with repeated train/test splits and 5-fold cross-validation; reports Test/CV R² and
MSE distributions, averaged feature importances, a feature-ablation analysis, and a
feature-set comparison. Input: `data/final_filtered_gene_tf_data.csv`.

---

## Reproducibility notes
- **Notebook outputs were stripped** and all absolute paths rewritten to be
  relative to `data/`, so notebooks are clean and portable.
- **18 referenced inputs were not located** in the source tree (GC/TATA Venn inputs
  and some RSEM DE tables) — listed in `data/README.md`. Regenerate them from the
  upstream step to reproduce those specific panels.
- **R package versions are unpinned** — after building the env, capture exact
  versions with `sessionInfo()` (see `environment/R_sessionInfo_notes.md`).

## Code availability
See [`CODE_AVAILABILITY.md`](CODE_AVAILABILITY.md) for a manuscript-ready statement.

## License & citation
MIT (see [`LICENSE`](LICENSE)). Please cite the associated manuscript and this
repository (see [`CITATION.cff`](CITATION.cff)).
