# Data

Processed input files consumed by the analysis code. **All code is designed to be
run with *this* `data/` directory as the working directory** — file paths in the
notebooks and R scripts are relative to here.

```bash
# Notebooks: launch Jupyter from data/ so the kernel's working dir is correct
cd data
jupyter lab           # then open ../bioinformatics/... or ../machine_learning/...

# R scripts: run them from data/ too
cd data
Rscript ../bioinformatics/plots_and_figures/mef_plots.R
```

`_MANIFEST.txt` lists every bundled file.

## Contents (by analysis stage)
- **DE inputs:** raw/normalized count matrices (`all_unnorm_counts_*`,
  `*_Deseq2_*normalized.csv`, `normalised_*_countdata.csv`), sample sheets
  (`MEF_counts/coldata_*`, `coldata_minus_het_lmc_*`), gene lists.
- **Motif inputs:** HOMER `knownResults_*` tables under
  `{MEF,LMC}_{I,II,III}_{up,down}regulated/`, `*_top20Homer_knownmotifs.csv`,
  TF-family motif summaries.
- **Enrichment inputs:** `*_GO_*_GSEA_clusterprofiler.csv`.
- **ML inputs:** `final_filtered_gene_tf_data.csv`, `merged_tf_family_counts.csv`,
  DE `*_rf.csv` lists, and `feb11_analysis/` (TF annotations, ENCODE `*.bed.gz`,
  `fimo_out/fimo.tsv`, promoter beds, `refGene.txt`).

## Missing inputs (18 files not located in the source tree)
These are referenced by the code but were not found when assembling the package;
regenerate them from the corresponding upstream step, or add them here before
reproducing the affected figures:

- GC/TATA-box Venn inputs: `GCbox_{MEF,LMC}_Dec6_2023.csv`,
  `TATA_{MEF,LMC}_Dec6_2023.csv`, `otherTFs_{MEF,LMC}_Dec6_2023.csv`,
  `outlier_genecard_genes_with_*_box_*_Dec{5,7}_2023.csv`
  → produced by the TATA/GC-box annotation step feeding
  `bioinformatics/plots_and_figures/GC_TATA_venndiagrams_MEF_LMC.ipynb`.
- RSEM DE result tables: `results_DEgenes_*_rsem_Dec8_2023.csv`
  → produced by the DESeq2 DE notebooks.

Also note: `bioinformatics/motif_analysis/TF_binding_sites_code.R` originally read
genecards inputs from `LMC_results/genecards/` and `MEF_results/`; place those
folders here under `data/` to reproduce that section.
