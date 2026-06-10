# Supplementary tables

Differential-expression and pathway-enrichment results for all six contrasts
(**LMC** and **MEF**, comparisons **I, II, III**). Provided both as per-comparison
CSVs and as two multi-sheet Excel workbooks (one sheet per comparison).

```
supplementary_tables/
├── Supplementary_Table_S1_DEgenes.xlsx        # 6 sheets: LMC_I/II/III, MEF_I/II/III
├── Supplementary_Table_S2_GSEA_pathways.xlsx  # 6 sheets: LMC_I/II/III, MEF_I/II/III
├── differential_expression/                   # per-comparison DE CSVs
└── gsea_pathway/                              # per-comparison GSEA CSVs
```

## Supplementary Table S1 — Differentially expressed genes

DESeq2 results, **significant genes only (adjusted p-value ≤ 0.05)**, sorted by
adjusted p-value. Columns:

| Column | Description |
|---|---|
| `gene` | Gene symbol |
| `log2FC` | log2 fold-change (Sp3-null vs. control) |
| `pvalue` | Wald test p-value |
| `adjusted_pvalue` | Benjamini–Hochberg adjusted p-value (DESeq2 `padj`) |

| Comparison | DE genes |
|---|---|
| LMC_I | 1,036 |
| LMC_II | 2,369 |
| LMC_III | 4,051 |
| MEF_I | 3,897 |
| MEF_II | 908 |
| MEF_III | 7,726 |

## Supplementary Table S2 — GSEA pathway enrichment (clusterProfiler)

Gene-set enrichment analysis over Gene Ontology terms (clusterProfiler `gseGO`),
run separately for each comparison, sorted by adjusted p-value. Columns:

| Column | Description |
|---|---|
| `ID` | GO term identifier |
| `Description` | GO term name |
| `setSize` | Number of genes in the term |
| `enrichmentScore` | Running enrichment score |
| `NES` | Normalized enrichment score |
| `pvalue` | GSEA p-value |
| `adjusted_pvalue` | Benjamini–Hochberg adjusted p-value (`p.adjust`) |
| `qvalue` | q-value |
| `leading_edge` | Leading-edge statistics |
| `core_enrichment` | Leading-edge (core) genes driving enrichment |

| Comparison | GO terms |
|---|---|
| LMC_I | 247 |
| LMC_II | 642 |
| LMC_III | 553 |
| MEF_I | 1,023 |
| MEF_II | 378 |
| MEF_III | 1,694 |

---
_Source: DESeq2 result tables (`results_DEgenes_*_2020counts*.csv`) and
clusterProfiler GSEA outputs (`{lmc,mef}_GO_{I,II,III}_GSEA_clusterprofiler.csv`)
from the analysis in `bioinformatics/`._
