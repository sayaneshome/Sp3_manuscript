setwd('Desktop/Bioinformatics_projects/Bulk_RNA_seqdata_sp3/cluster_data/')

library(pheatmap)

csv1 <- read.csv('filter_MEF.csv')
colnames(csv1)[1] <- "Gene_name"
rownames(csv1) <- csv1$Gene_name
csv1$Gene_name <- NULL
csv1_m <- as.matrix(csv1)


annotation_colors <- list(
Genotype = c(WT = "#0bd8e0", KO = "#d8a2eb")  # Colors for both MEF and LMC
)

genotype_info_MEF <- data.frame(
Genotype = c("WT", "WT","WT", "WT", "WT","KO", "KO", "KO")
)
rownames(genotype_info_MEF) <- colnames(csv1_m)

heatmap_result_csv1 <- pheatmap(csv1_m,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 0.75,
         annotation_col = genotype_info_MEF,
         annotation_colors = annotation_colors)
# Extract gene order from MEF heatmap
gene_order_MEF <- rownames(csv1_m)[heatmap_result_csv1$tree_row$order]

csv2 <- read.csv('LMC_Deseq2_countsdata_normalized.csv')
cols_to_remove <- grep("LPS", names(csv2), ignore.case = TRUE)
csv2 <- csv2[,-cols_to_remove]
rownames(csv2) <- csv2$X
csv2$X <- NULL

csv2_filtered <- csv2[rownames(csv1), ]  # Ensure LMC data includes genes from MEF data
csv2_ordered <- csv2_filtered[gene_order_MEF, ]
csv2_cleaned <- na.omit(csv2_ordered)

# Ensure that only row names present in both csv1_m and csv2_cleaned are used
row_names_to_keep <- intersect(rownames(csv1_m), rownames(csv2_cleaned))

# Now, filter csv1_m to keep only those rows that exist in csv2_cleaned's row names
csv1_m_filtered <- csv1_m[row_names_to_keep, , drop = FALSE]
csv2_m_filtered <- csv2_cleaned[row_names_to_keep, , drop = FALSE]

chk <- pheatmap(csv1_m_filtered,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 0.25,
         annotation_col = genotype_info_MEF,
         annotation_colors = annotation_colors)

gene_order_MEF <- rownames(csv1_m_filtered)[chk$tree_row$order]

genotype_info_LMC <- data.frame(
  Genotype = c("WT", "WT","WT", "WT", "KO","KO", "KO", "KO")
)
rownames(genotype_info_LMC) <- colnames(csv2_m_filtered)

csv2_m_filtered <- as.matrix(csv2_m_filtered)#have to work on this

# Assuming gene_order_MEF is already defined and contains the desired gene order
# Ensure csv2_m_filtered has rows in the same order as gene_order_MEF
csv2_m_ordered <- csv2_m_filtered[match(gene_order_MEF, rownames(csv2_m_filtered)),]


pheatmap(csv2_m_ordered,
scale = "row",
cluster_rows = FALSE,
cluster_cols = FALSE,
show_rownames = TRUE,
show_colnames = FALSE,
fontsize_row = 0.25,
annotation_col = genotype_info_LMC,
annotation_colors = annotation_colors)

