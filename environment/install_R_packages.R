# R package installation for the Sp3 RNA-seq reproducibility package
# Developed against R 4.4.x.  Run:  Rscript install_R_packages.R
# These packages are loaded by the differential_expression notebooks (via rpy2 %%R)
# and by the standalone .R scripts in bioinformatics/.

# ---- CRAN packages ----
cran <- c(
  "tidyverse", "dplyr", "tidyr", "stringr", "tibble", "readxl", "reshape2",
  "ggplot2", "ggrepel", "ggridges", "cowplot", "gridExtra", "patchwork",
  "pheatmap", "RColorBrewer", "corrplot", "factoextra", "VennDiagram",
  "UpSetR", "ComplexUpset", "plotly", "svglite", "systemfonts",
  "randomForest", "ashr", "jsonlite", "htmltools", "knitr", "repr"
)
new_cran <- cran[!cran %in% rownames(installed.packages())]
if (length(new_cran)) install.packages(new_cran, repos = "https://cloud.r-project.org")

# ---- Bioconductor packages ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

bioc <- c(
  "DESeq2", "apeglm",
  "clusterProfiler", "enrichplot", "EnhancedVolcano", "Glimma",
  "org.Mm.eg.db", "org.Hs.eg.db", "EnsDb.Mmusculus.v79"
)
new_bioc <- bioc[!bioc %in% rownames(installed.packages())]
if (length(new_bioc)) BiocManager::install(new_bioc, update = FALSE, ask = FALSE)

cat("\nDone. Run sessionInfo() after loading the libraries to record exact versions.\n")
