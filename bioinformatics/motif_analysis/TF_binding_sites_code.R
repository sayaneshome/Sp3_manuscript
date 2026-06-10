# Run this script from the repository's data/ directory as the working directory.
# (Originally pointed at .../cluster_data/LMC_results/genecards/ — supply those
#  genecards inputs under data/LMC_results/genecards/ if reproducing this section.)
# setwd('LMC_results/genecards/')
library(readxl)
library(dplyr)

#analysing TF-binding sites from DE genes from LMC dataset
file1 <- read_xlsx('genecards_file1.xlsx',4)
file2 <- read_xlsx('genecards_file2.xlsx',3)
file3 <- read_xlsx('genecards_file3.xlsx',5)
file4 <- read_xlsx('genecards_file4.xlsx',4)
file5 <- read_xlsx('genecards_file5.xlsx',4)
file6 <- read_xlsx('genecards_file6.xlsx',6)
file7 <- read_xlsx('genecards_file7.xlsx',3)
file1_short <- data.frame(file1$Symbol,file1$TFBSs)
file2_short <- data.frame(file2$Symbol,file2$TFBSs)
file3_short <- data.frame(file3$Symbol,file3$TFBSs)
file4_short <- data.frame(file4$Symbol,file4$TFBSs)
file5_short <- data.frame(file5$Symbol,file5$TFBSs)
file6_short <- data.frame(file6$Symbol,file6$TFBSs)
file7_short <- data.frame(file7$Symbol,file7$TFBSs)
colnames(file1_short)[1] <- "Gene_name"
colnames(file1_short)[2] <- "TFBS"
colnames(file2_short)[1] <- "Gene_name"
colnames(file3_short)[1] <- "Gene_name"
colnames(file4_short)[1] <- "Gene_name"
colnames(file5_short)[1] <- "Gene_name"
colnames(file6_short)[1] <- "Gene_name"
colnames(file7_short)[1] <- "Gene_name"
colnames(file2_short)[2] <- "TFBS"
colnames(file3_short)[2] <- "TFBS"
colnames(file4_short)[2] <- "TFBS"
colnames(file5_short)[2] <- "TFBS"
colnames(file6_short)[2] <- "TFBS"
colnames(file7_short)[2] <- "TFBS"

TF_df <- rbind(file1_short,file2_short,file3_short,file4_short,file5_short,file6_short,file7_short)

library(tidyverse)

# Separate values into different rows and remove NAs
separated_df <- TF_df %>%
  separate_rows(TFBS, sep = "\\|\\|") %>%
  na.omit()

# Print the resulting dataframe
znf_rows <- separated_df[grepl("^ZNF|^SP|^KLF", separated_df$TFBS), ]
znf_rows <- unique(znf_rows)
tbp_rows <- subset(separated_df, TFBS == "TBP")
tbp_rows <- unique(tbp_rows)

promoters <- separated_df[!grepl("TBP|ZNF|SP|KLF", separated_df$TFBS), ]
promoters <- unique(promoters)

write.csv(tbp_rows,'tata_box_binding_sites_lmc.csv')
write.csv(znf_rows,'gc_box_binding_sites_lmc.csv')
write.csv(promoters,'TF_binding_sites_lmc.csv')

# setwd('MEF_results/')   # supply MEF genecards inputs under data/MEF_results/ to reproduce this section
library(readxl)
mef_1 <- read_xlsx('results_DEgenes_mut_LPStreatment_MEF.xlsx')
mef_2 <- read_xlsx('results_DEgenes_mutant_versus_wildtype_considering_LPStreatment_MEF.xlsx')
mef_3 <- read_xlsx('results_DEgenes_wildtype_sp3knockout_onLPStreatment_MEF.xlsx')
mef_4 <- read_xlsx('results_DEgenes_wt_LPStreatment_MEF.xlsx')
mef_5 <- read_xlsx('results_DEgenes_WT_Sp3knockout_global_MEF.xlsx')

mef_1 <- mef_1[3:nrow(mef_1),]
mef_2 <- mef_2[3:nrow(mef_2),]
mef_3 <- mef_3[2:nrow(mef_3),]
mef_4 <- mef_4[2:nrow(mef_4),]
mef_5 <- mef_5[3:nrow(mef_4),]

colnames(mef_1) <- as.character(mef_1[1, ])
# Remove the first row from the dataframe
mef_1 <- mef_1[-1, ]
colnames(mef_2) <- as.character(mef_2[1, ])
# Remove the first row from the dataframe
mef_2 <- mef_2[-1, ]
colnames(mef_3) <- as.character(mef_3[1, ])
# Remove the first row from the dataframe
mef_3 <- mef_3[-1, ]
colnames(mef_4) <- as.character(mef_4[1, ])
# Remove the first row from the dataframe
mef_4 <- mef_4[-1, ]
colnames(mef_5) <- as.character(mef_5[1, ])
# Remove the first row from the dataframe
mef_5 <- mef_5[-1, ]
mef_1_genes <- data.frame(mef_1$`Gene name`)
mef_2_genes <- data.frame(mef_2$`Gene name`)
mef_3_genes <- data.frame(mef_3$SYMBOL)
mef_4_genes <- data.frame(mef_4$`Gene name`)
mef_5_genes <- data.frame(mef_5$`Gene name`)
colnames(mef_1_genes)[1] <- "Gene_name"
colnames(mef_2_genes)[1] <- "Gene_name"
colnames(mef_3_genes)[1] <- "Gene_name"
colnames(mef_4_genes)[1] <- "Gene_name"
colnames(mef_5_genes)[1] <- "Gene_name"
mef_genes <- rbind(mef_1_genes,mef_2_genes,mef_3_genes,mef_4_genes,mef_5_genes)
mef_genes$Gene_name <- toupper(mef_genes$Gene_name)
write.csv(mef_genes,'mef_genes.csv')
mef_genes1 <- read.csv('mef_genes.csv',header = FALSE)
colnames(mef_genes1)[1] <- "Gene_name"
tf_mef <- merge(mef_genes1,separated_df,by = "Gene_name")
tf_mef <- unique(tf_mef)

mef_genes1$Gene_name <- as.character(mef_genes1$Gene_name)
separated_df$Gene_name <- as.character(separated_df$Gene_name)
result <- anti_join(mef_genes1, separated_df, by = c("Gene_name"))
genenames_lmc <- read.csv('LMC_results/gene_names.csv')
genenames_lmc$X <- NULL
colnames(genenames_lmc)[1] <- "Gene_name"
genenames_mef <- read.csv('mef_genes.csv',header = FALSE)
genenames_lmc <- read.csv('LMC_results/gene_names.csv',header = FALSE)
colnames(genenames_lmc)[1] <- "Gene_name"
colnames(genenames_mef)[1] <- "Gene_name"
genenames_lmc$Gene_name <- toupper(genenames_lmc$Gene_name)
genenames_mef$Gene_name <- as.character(genenames_mef$Gene_name)
genenames_lmc$Gene_name <- as.character(genenames_lmc$Gene_name)

result_mef <- anti_join(genenames_mef, genenames_lmc, by = "Gene_name")
result_lmc <- anti_join(genenames_lmc, genenames_mef, by = "Gene_name")
write.csv(result_mef,"DE_genes_in_MEFdataset_not_in_LMCdataset.csv")
write.csv(result_lmc,"DE_genes_in_LMCdataset_not_in_MEFdataset.csv")
filtered_df_mef <- tf_mef[!grepl("TBP|ZNF", tf_mef$TFBS), ]
# Print the filtered dataframe
tbp_mef <- subset(tf_mef, TFBS == "TBP")
gc_box_mef <- tf_mef[grepl("^ZNF|^SP|^KLF", tf_mef$TFBS), ]
common_names <- merge(genenames_lmc,genenames_mef,by = "Gene_name")
write.csv(common_names,"DEgenes_in_both_LMC_and_MEF_datasets.csv")

#Updating MEF list
csv1 <- read_xlsx('GeneALaCart-514864-230526-124453 (1).xlsx',3)
csv2 <- read_xlsx('GeneALaCart-514864-230526-182934.xlsx',3)

csv1_short <- data.frame(csv1$Symbol,csv1$TFBSs)
csv2_short <- data.frame(csv2$Symbol,csv2$TFBSs)
colnames(csv1_short)[1] <- "Gene_name"
colnames(csv1_short)[2] <- "TFBS"
colnames(csv2_short)[1] <- "Gene_name"
colnames(csv2_short)[2] <- "TFBS"
Tf_df1 <- rbind(csv1_short,csv2_short)
# Separate values into different rows and remove NAs
separated_df1 <- Tf_df1 %>%
  separate_rows(TFBS, sep = "\\|\\|") %>%
  na.omit()


znf_rows <- separated_df1[grepl("^ZNF|^KLF|^SP", separated_df$TFBS), ]
znf_rows <- unique(znf_rows)
gc_box_mef <- unique(gc_box_mef)
gc_box_mef1 <- rbind(gc_box_mef,znf_rows)


tbp_rows <- subset(separated_df1, TFBS == "TBP")
tbp_rows <- unique(tbp_rows$Gene_name)
tbp_rows <- subset(separated_df1, TFBS == "TBP")
tbp_rows <- unique(tbp_rows)

znf_rows <- unique(znf_rows)
tbp_mef <- unique(tbp_mef)
tbp_mef1 <- rbind(znf_rows,tbp_mef)

promoters <- separated_df1[!grepl("TBP|ZNF|KLF|SP", separated_df1$TFBS), ]
promoters <- unique(promoters)
filtered_df_mef <- unique(filtered_df_mef)
filtered_df_mef1 <- rbind(filtered_df_mef,promoters)

write.csv(gc_box_mef1,"gc_box_binding_sites_mef.csv")
write.csv(filtered_df_mef1,'TF_binding_sites_mef.csv')
write.csv(tbp_mef1,"tata_box_binding_sites_mef.csv")