# R environment notes

The R analyses (DESeq2 differential expression, clusterProfiler GO/GSEA,
EnhancedVolcano, Glimma, motif/heatmap plotting) were developed on **R 4.4.x**.

Exact package versions were **not** pinned in the original environment, so they
cannot be reconstructed from the code alone. To record the versions you actually
use for the submission, run this once after installing (see `install_R_packages.R`)
and loading the libraries in a notebook/session:

```r
sink("environment/R_sessionInfo.txt")
sessionInfo()
sink()
```

Commit the resulting `R_sessionInfo.txt` alongside this file so reviewers have the
exact build.

## Packages used (loaded across the notebooks and .R scripts)

**Differential expression / annotation (Bioconductor):**
DESeq2, apeglm, ashr, EnhancedVolcano, Glimma, clusterProfiler, enrichplot,
org.Mm.eg.db, org.Hs.eg.db, EnsDb.Mmusculus.v79

**Plotting / tables (CRAN):**
tidyverse, dplyr, tidyr, stringr, tibble, readxl, reshape2, ggplot2, ggrepel,
ggridges, cowplot, gridExtra, patchwork, pheatmap, RColorBrewer, corrplot,
factoextra, VennDiagram, UpSetR, ComplexUpset, plotly, svglite, systemfonts,
randomForest, jsonlite, htmltools, knitr, repr
