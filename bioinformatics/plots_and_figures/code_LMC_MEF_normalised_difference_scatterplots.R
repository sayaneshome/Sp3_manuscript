setwd('/Users/sshome/Desktop/Bioinformatics_projects/Bulk_RNA_seqdata_sp3/cluster_data/')
LMC_norm <- read.csv('LMC_normalized_counts_differences.csv')
MEF_norm <- read.csv('MEF_normalized_counts_differences.csv')
LMC_norm$mean_mut <- (LMC_norm$X697.5 + LMC_norm$X697.2 + LMC_norm$X685.1 + LMC_norm$X468.6)/4
View(LMC_norm)
LMC_norm$mean_wt <- (LMC_norm$X468.4 + LMC_norm$X685.2 + LMC_norm$X697.1 + LMC_norm$X697.3)/4
LMC_norm$difference_mean <- - LMC_norm$mean_mut + LMC_norm$mean_wt
LMC_norm$X <- NULL
View(LMC_norm)

MEF_norm$mean_mut <- (MEF_norm$P4_5_1 + MEF_norm$P4_5_2 + MEF_norm$P4_5_3)/3
MEF_norm$mean_wt <- (MEF_norm$P7_4_1 + MEF_norm$P7_4_2 + MEF_norm$P8_4_3 + MEF_norm$P8_4_4 + MEF_norm$P8_4_5)/5
MEF_norm$difference_mean <- MEF_norm$mean_wt - MEF_norm$mean_mut  
MEF_norm$X <- NULL


setwd('/Users/sshome/Desktop/Bioinformatics_projects/Bulk_RNA_seqdata_sp3/cluster_data/Homer_results/Homer_motifs_known_shortlisted/')

homer_I_MEF_down <- read.csv('mef_I_downregulated_top20Homer_knownmotifs.csv')
homer_I_MEF_up <- read.csv('mef_I_upregulated_top20Homer_knownmotifs.csv')
homer_II_MEF_up <- read.csv('mef_II_upregulated_top20Homer_knownmotifs.csv')
homer_II_MEF_down <- read.csv('mef_II_downregulated_top20Homer_knownmotifs.csv')
homer_III_MEF_down <- read.csv('mef_III_downregulated_top20Homer_knownmotifs.csv')
homer_III_MEF_up <- read.csv('mef_III_upregulated_top20Homer_knownmotifs.csv')
homer_I_LMC_down <- read.csv('lmc_I_downregulated_top20Homer_knownmotifs.csv')
homer_I_LMC_up <- read.csv('lmc_I_upregulated_top20Homer_knownmotifs.csv')
homer_II_LMC_up <- read.csv('lmc_II_upregulated_top20Homer_knownmotifs.csv')
homer_II_LMC_down <- read.csv('lmc_II_downregulated_top20Homer_knownmotifs.csv')
homer_III_LMC_up <- read.csv('lmc_III_upregulated_top20Homer_knownmotifs.csv')
homer_III_LMC_down <- read.csv('lmc_III_downregulated_top20Homer_knownmotifs.csv')


homer_I_LMC_down$X <- NULL
homer_II_LMC_down$X <- NULL
homer_III_LMC_down$X <- NULL
homer_III_LMC_up$X <- NULL
homer_II_LMC_up$X <- NULL
homer_I_LMC_up$X <- NULL
homer_I_MEF_up$X <- NULL
homer_II_MEF_up$X <- NULL
homer_III_MEF_up$X <- NULL
homer_III_MEF_down$X <- NULL
homer_II_MEF_down$X <- NULL
homer_I_MEF_down$X <- NULL

library(dplyr)

homer_I_LMC_down<- homer_I_LMC_down %>% filter(Motif.Name != "unknown")
homer_I_LMC_up<- homer_I_LMC_up %>% filter(Motif.Name != "unknown")
homer_II_LMC_up<- homer_II_LMC_up %>% filter(Motif.Name != "unknown")
homer_II_LMC_down<- homer_II_LMC_down %>% filter(Motif.Name != "unknown")
homer_III_LMC_down<- homer_III_LMC_down %>% filter(Motif.Name != "unknown")
homer_III_LMC_up<- homer_III_LMC_up %>% filter(Motif.Name != "unknown")
homer_III_MEF_up<- homer_III_MEF_up %>% filter(Motif.Name != "unknown")
homer_III_MEF_down<- homer_III_MEF_down %>% filter(Motif.Name != "unknown")
homer_II_MEF_down<- homer_II_MEF_down %>% filter(Motif.Name != "unknown")
homer_II_MEF_up<- homer_II_MEF_up %>% filter(Motif.Name != "unknown")
homer_I_MEF_up<- homer_I_MEF_up %>% filter(Motif.Name != "unknown")
homer_I_MEF_down<- homer_I_MEF_down %>% filter(Motif.Name != "unknown")


homer_I_LMC_down <- homer_I_LMC_down %>% select(Motif.Name, Consensus,Log.P.value,q.value..Benjamini.,X..of.Target.Sequences.with.Motif.1,X..of.Background.Sequences.with.Motif.1)
homer_II_LMC_down <- homer_II_LMC_down %>% select(Motif.Name, Consensus,Log.P.value,q.value..Benjamini.,X..of.Target.Sequences.with.Motif.1,X..of.Background.Sequences.with.Motif.1)
homer_III_LMC_down <- homer_III_LMC_down %>% select(Motif.Name, Consensus,Log.P.value,q.value..Benjamini.,X..of.Target.Sequences.with.Motif.1,X..of.Background.Sequences.with.Motif.1)
homer_III_LMC_up <- homer_III_LMC_up %>% select(Motif.Name, Consensus,Log.P.value,q.value..Benjamini.,X..of.Target.Sequences.with.Motif.1,X..of.Background.Sequences.with.Motif.1)
homer_II_LMC_up <- homer_II_LMC_up %>% select(Motif.Name, Consensus,Log.P.value,q.value..Benjamini.,X..of.Target.Sequences.with.Motif.1,X..of.Background.Sequences.with.Motif.1)
homer_I_LMC_up <- homer_I_LMC_up %>% select(Motif.Name, Consensus,Log.P.value,q.value..Benjamini.,X..of.Target.Sequences.with.Motif.1,X..of.Background.Sequences.with.Motif.1)
homer_I_MEF_up <- homer_I_MEF_up %>% select(Motif.Name, Consensus,Log.P.value,q.value..Benjamini.,X..of.Target.Sequences.with.Motif.1,X..of.Background.Sequences.with.Motif.1)
homer_II_MEF_up <- homer_II_MEF_up %>% select(Motif.Name, Consensus,Log.P.value,q.value..Benjamini.,X..of.Target.Sequences.with.Motif.1,X..of.Background.Sequences.with.Motif.1)
homer_III_MEF_up <- homer_III_MEF_up %>% select(Motif.Name, Consensus,Log.P.value,q.value..Benjamini.,X..of.Target.Sequences.with.Motif.1,X..of.Background.Sequences.with.Motif.1)
homer_III_MEF_down <- homer_III_MEF_down %>% select(Motif.Name, Consensus,Log.P.value,q.value..Benjamini.,X..of.Target.Sequences.with.Motif.1,X..of.Background.Sequences.with.Motif.1)
homer_II_MEF_down <- homer_II_MEF_down %>% select(Motif.Name, Consensus,Log.P.value,q.value..Benjamini.,X..of.Target.Sequences.with.Motif.1,X..of.Background.Sequences.with.Motif.1)
homer_I_MEF_down <- homer_I_MEF_down %>% select(Motif.Name, Consensus,Log.P.value,q.value..Benjamini.,X..of.Target.Sequences.with.Motif.1,X..of.Background.Sequences.with.Motif.1)

# List all objects in the environment that start with "Homer"
homer_objects <- ls(pattern = "^homer")

# Get the actual objects
homer_dataframes <- mget(homer_objects)

# Define new column names
new_col_names <- c("Motif Name", "Consensus", "Log P-value", 
                   "q-value (Benjamini)", "% of Target Sequences with Motif", 
                   "% of Background Sequences with Motif")

# Loop through each dataframe and update column names
for (df_name in names(homer_dataframes)) {
  current_df <- homer_dataframes[[df_name]]
  
  # Check if it's indeed a dataframe and has at least 6 columns
  if (is.data.frame(current_df) && ncol(current_df) >= 6) {
    colnames(current_df)[1:6] <- new_col_names
    
    # Assign the modified dataframe back to the environment
    assign(df_name, current_df, envir = .GlobalEnv)
  }
}

# Initialize an empty list to store the filtered dataframes
filtered_dfs <- list()
# List all objects in the environment that start with "Homer"
homer_objects <- ls(pattern = "^homer")
# Loop through each object
for (df_name in homer_objects) {
# Use get() to retrieve the object based on its name
current_df <- get(df_name)
# Ensure it's a dataframe and has the 'Motif Name' column
if (is.data.frame(current_df) && "Motif Name" %in% colnames(current_df)) {
# Filter rows where 'Motif Name' contains 'ZF'
matching_rows <- filter(current_df, grepl("Zf", `Motif Name`))
# Check if any rows match the criteria
if (nrow(matching_rows) > 0) {
# Add a new column to indicate the source DataFrame
matching_rows$Source_DF <- df_name
# Add the modified dataframe to the list
filtered_dfs[[df_name]] <- matching_rows
}
}
}
# Combine all filtered DataFrames into one
if (length(filtered_dfs) > 0) {
combined_df <- do.call("rbind", filtered_dfs)
# Reset row names
row.names(combined_df) <- NULL
} else {
combined_df <- NULL
print("No rows containing 'ZF' were found in 'Motif Name' columns.")
}
# Check the final combined dataframe
if (!is.null(combined_df)) {
print("Combined dataframe created successfully.")
} else {
print("No matching data found across all DataFrames.")
}
View(combined_df)
ZF_motifs_combined <- combined_df
# Initialize an empty list to store the filtered dataframes
filtered_dfs <- list()
# List all objects in the environment that start with "Homer"
homer_objects <- ls(pattern = "^homer")
# Loop through each object
for (df_name in homer_objects) {
# Use get() to retrieve the object based on its name
current_df <- get(df_name)
# Ensure it's a dataframe and has the 'Motif Name' column
if (is.data.frame(current_df) && "Motif Name" %in% colnames(current_df)) {
# Filter rows where 'Motif Name' contains 'Homeobox'
matching_rows <- filter(current_df, grepl("Homeobox", `Motif Name`))
# Check if any rows match the criteria
if (nrow(matching_rows) > 0) {
# Add a new column to indicate the source DataFrame
matching_rows$Source_DF <- df_name
# Add the modified dataframe to the list
filtered_dfs[[df_name]] <- matching_rows
}
}
}
# Combine all filtered DataFrames into one
if (length(filtered_dfs) > 0) {
combined_df <- do.call("rbind", filtered_dfs)
# Reset row names
row.names(combined_df) <- NULL
} else {
combined_df <- NULL
print("No rows containing 'Homeobox' were found in 'Motif Name' columns.")
}
# Check the final combined dataframe
if (!is.null(combined_df)) {
print("Combined dataframe created successfully.")
} else {
print("No matching data found across all DataFrames.")
}
View(combined_df)
Homeobox_motifs_combined <- combined_df
View(homer_I_LMC_down)
# Initialize an empty list to store the filtered dataframes
filtered_dfs <- list()
# List all objects in the environment that start with "Homer"
homer_objects <- ls(pattern = "^homer")
# Loop through each object
for (df_name in homer_objects) {
# Use get() to retrieve the object based on its name
current_df <- get(df_name)
# Ensure it's a dataframe and has the 'Motif Name' column
if (is.data.frame(current_df) && "Motif Name" %in% colnames(current_df)) {
# Filter rows where 'Motif Name' contains 'IRF'
matching_rows <- filter(current_df, grepl("IRF", `Motif Name`))
# Check if any rows match the criteria
if (nrow(matching_rows) > 0) {
# Add a new column to indicate the source DataFrame
matching_rows$Source_DF <- df_name
# Add the modified dataframe to the list
filtered_dfs[[df_name]] <- matching_rows
}
}
}
# Combine all filtered DataFrames into one
if (length(filtered_dfs) > 0) {
combined_df <- do.call("rbind", filtered_dfs)
# Reset row names
row.names(combined_df) <- NULL
} else {
combined_df <- NULL
print("No rows containing 'IRF' were found in 'Motif Name' columns.")
}
# Check the final combined dataframe
if (!is.null(combined_df)) {
print("Combined dataframe created successfully.")
} else {
print("No matching data found across all DataFrames.")
}
View(combined_df)
IRF_motifs_combined <- combined_df
View(homer_I_LMC_down)
# Initialize an empty list to store the filtered dataframes
filtered_dfs <- list()
# List all objects in the environment that start with "Homer"
homer_objects <- ls(pattern = "^homer")
# Loop through each object
for (df_name in homer_objects) {
# Use get() to retrieve the object based on its name
current_df <- get(df_name)
# Ensure it's a dataframe and has the 'Motif Name' column
if (is.data.frame(current_df) && "Motif Name" %in% colnames(current_df)) {
# Filter rows where 'Motif Name' contains 'NR'
matching_rows <- filter(current_df, grepl("NR", `Motif Name`))
# Check if any rows match the criteria
if (nrow(matching_rows) > 0) {
# Add a new column to indicate the source DataFrame
matching_rows$Source_DF <- df_name
# Add the modified dataframe to the list
filtered_dfs[[df_name]] <- matching_rows
}
}
}
# Combine all filtered DataFrames into one
if (length(filtered_dfs) > 0) {
combined_df <- do.call("rbind", filtered_dfs)
# Reset row names
row.names(combined_df) <- NULL
} else {
combined_df <- NULL
print("No rows containing 'NR' were found in 'Motif Name' columns.")
}
# Check the final combined dataframe
if (!is.null(combined_df)) {
print("Combined dataframe created successfully.")
} else {
print("No matching data found across all DataFrames.")
}
View(combined_df)
NR_motifs_combined <- combined_df
View(homer_I_LMC_down)
# Initialize an empty list to store the filtered dataframes
filtered_dfs <- list()
# List all objects in the environment that start with "Homer"
homer_objects <- ls(pattern = "^homer")
# Loop through each object
for (df_name in homer_objects) {
# Use get() to retrieve the object based on its name
current_df <- get(df_name)
# Ensure it's a dataframe and has the 'Motif Name' column
if (is.data.frame(current_df) && "Motif Name" %in% colnames(current_df)) {
# Filter rows where 'Motif Name' contains 'Stat'
matching_rows <- filter(current_df, grepl("Stat", `Motif Name`))
# Check if any rows match the criteria
if (nrow(matching_rows) > 0) {
# Add a new column to indicate the source DataFrame
matching_rows$Source_DF <- df_name
# Add the modified dataframe to the list
filtered_dfs[[df_name]] <- matching_rows
}
}
}
# Combine all filtered DataFrames into one
if (length(filtered_dfs) > 0) {
combined_df <- do.call("rbind", filtered_dfs)
# Reset row names
row.names(combined_df) <- NULL
} else {
combined_df <- NULL
print("No rows containing 'Stat' were found in 'Motif Name' columns.")
}
# Check the final combined dataframe
if (!is.null(combined_df)) {
print("Combined dataframe created successfully.")
} else {
print("No matching data found across all DataFrames.")
}
STAT_motifs_combined <- combined_df

# List all objects in the environment that end with 'combined'
combined_objects <- ls(pattern = "combined$")

# Loop through each object name
for (df_name in combined_objects) {
  # Retrieve the object based on its name
  current_df <- get(df_name)
  
  # Ensure it's a dataframe
  if (is.data.frame(current_df)) {
    # Ensure the dataframe has at least 6 columns
    if (ncol(current_df) >= 6) {
      # Convert 5th and 6th columns by removing "%" and converting to numeric
      current_df[,5] <- as.numeric(gsub("%", "", current_df[,5]))
      current_df[,6] <- as.numeric(gsub("%", "", current_df[,6]))
      
      # Sort the dataframe by the 5th column in descending order
      current_df <- current_df[order(current_df[,5], decreasing = TRUE), ]
      
      # Assign the modified and sorted dataframe back to the global environment with the original name
      assign(df_name, current_df, envir = .GlobalEnv)
    } else {
      message(paste(df_name, "does not have at least 6 columns. Skipping."))
    }
  } else {
    message(paste(df_name, "is not a dataframe. Skipping."))
  }
}

print("Dataframes ending with 'combined' have been updated: 5th and 6th columns had '%' removed, converted to numeric, then sorted by the 5th column in descending order.")



write.csv(IRF_motifs_combined,"all_IRF_family_motifs_combined_homeranalysis.csv")
write.csv(ZF_motifs_combined,"all_ZF_family_motifs_combined_homeranalysis.csv")
write.csv(NR_motifs_combined,"all_NR_family_motifs_combined_homeranalysis.csv")
write.csv(STAT_motifs_combined,"all_STAT_family_motifs_combined_homeranalysis.csv")




