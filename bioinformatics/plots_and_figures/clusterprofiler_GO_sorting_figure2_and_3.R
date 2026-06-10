library(dplyr)
library(tidyr)
lmc_I_GO <- read.csv('lmc_GO_I_GSEA_clusterprofiler.csv')
lmc_II_GO <- read.csv('lmc_GO_II_GSEA_clusterprofiler.csv')

# Merge the dataframes
merged_df <- merge(lmc_I_GO, lmc_II_GO, by="ID", suffixes=c("_I", "_II"))

# Select, rename the relevant columns, and add the NES difference column
final_df <- merged_df %>%
  select(ID, Description_I, NES_I, NES_II) %>%
  rename(Description = Description_I,
         NES_LMC_I = NES_I,
         NES_LMC_II = NES_II) %>%
  mutate(NES_Difference = NES_LMC_I - NES_LMC_II)

# Sort by the absolute value of NES_Difference to highlight the largest differences
final_lmc_df <- final_df %>%
  arrange(desc(abs(NES_Difference)))

# Display the first few rows of the new dataframe
head(final_lmc_df, 10)

write.csv(final_lmc_df,"LMC_I_II_GO_comparison_figure2_3.csv")


library(dplyr)

mef_I_GO <- read.csv('mef_GO_I_GSEA_clusterprofiler.csv')
mef_II_GO <- read.csv('mef_GO_II_GSEA_clusterprofiler.csv')

# Merge the MEF dataframes
merged_mef_df <- merge(mef_I_GO, mef_II_GO, by="ID", suffixes=c("_I", "_II"))

# Select, rename the relevant columns, and add the NES difference column
final_mef_df <- merged_mef_df %>%
  select(ID, Description_I, NES_I, NES_II) %>%
  rename(Description = Description_I,
         NES_MEF_I = NES_I,
         NES_MEF_II = NES_II) %>%
  mutate(NES_Difference = NES_MEF_I - NES_MEF_II)

# Sort by the absolute value of NES_Difference to highlight the largest differences
final_mef_df <- final_mef_df %>%
  arrange(desc(abs(NES_Difference)))

# Display the first few rows of the new dataframe
head(final_mef_df, 10)

write.csv(final_mef_df,"MEF_I_II_GO_comparison_figure2_3.csv")
