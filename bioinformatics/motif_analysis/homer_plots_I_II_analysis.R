# Set working directory
setwd('Desktop/Bioinformatics_projects/Bulk_RNA_seqdata_sp3/cluster_data/Homer_results/Homer_motifs_known_shortlisted/')

# Load the necessary libraries
library(ggplot2)
library(gridExtra)
library(cowplot)
library(dplyr)
library(stringr)

# Define a mapping of old condition names to the new labels
condition_labels <- c(
  "I_MEF_up" = "Higher in wildtype [MEF]",
  "I_MEF_down" = "Lower in wildtype [MEF]",
  "I_LMC_up" = "Higher in wildtype [LMC]",
  "I_LMC_down" = "Lower in wildtype [LMC]",
  "II_MEF_up" = "Higher in mutant [MEF]",
  "II_MEF_down" = "Lower in mutant [MEF]",
  "II_LMC_up" = "Higher in mutant [LMC]",
  "II_LMC_down" = "Lower in mutant [LMC]"
)

# Function to extract the condition from the source dataframe strings
extract_condition <- function(df) {
  df$Condition <- str_extract(df$Source_DF, paste(names(condition_labels), collapse = "|"))
  df$Condition <- factor(df$Condition, levels = names(condition_labels))
  return(df)
}

# Function to filter dataframe for specified conditions
filter_conditions <- function(df, conditions) {
  df <- df %>% filter(Condition %in% conditions)
  return(df)
}

# Function to calculate frequency of each value in the condition column, ensuring all levels are included
calculate_frequency <- function(df) {
  freq_table <- table(df$Condition)
  freq_df <- as.data.frame(freq_table)
  colnames(freq_df) <- c("Condition", "Frequency")
  return(freq_df)
}

# Function to create bar plots
create_plot <- function(freq_df, title) {
  ggplot(freq_df, aes(x = Condition, y = Frequency, fill = Condition)) +
    geom_bar(stat = "identity", width = 0.5) +
    coord_flip() +
    ggtitle(title) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(hjust = 1, color = "black", size = 12),
      axis.text.y = element_text(color = "black", size = 10, margin = margin(r = 5)),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      legend.position = "none",
      plot.background = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)  # Reduce white space around the plot
    ) +
    geom_text(aes(label = Frequency), hjust = -0.2, size = 4) + # Add labels
    scale_x_discrete(labels = condition_labels[as.character(freq_df$Condition)]) +
    scale_fill_manual(values = c(
      "I_MEF_up" = "#96cb01", 
      "I_MEF_down" = "#fe9388", 
      "I_LMC_up" = "#96cb01", 
      "I_LMC_down" = "#fe9388",
      "II_MEF_up" = "#96cb01", 
      "II_MEF_down" = "#fe9388", 
      "II_LMC_up" = "#96cb01", 
      "II_LMC_down" = "#fe9388"
    )) +
    expand_limits(y = 0)  # Ensure bars start from the axis
}

# Read the CSV files
csv_homeobox <- read.csv('all_homeobox_family_motifs_combined_homeranalysis.csv')
csv_irf <- read.csv('all_IRF_family_motifs_combined_homeranalysis.csv')
csv_stat <- read.csv('all_STAT_family_motifs_combined_homeranalysis.csv')

# Extract conditions
csv_homeobox <- extract_condition(csv_homeobox)
csv_irf <- extract_condition(csv_irf)
csv_stat <- extract_condition(csv_stat)

# Filter dataframes for specified conditions (Analysis I and II)
conditions_to_plot <- c("I_MEF_up", "I_MEF_down", "I_LMC_up", "I_LMC_down", 
                        "II_MEF_up", "II_MEF_down", "II_LMC_up", "II_LMC_down")

csv_homeobox_filtered <- filter_conditions(csv_homeobox, conditions_to_plot)
csv_irf_filtered <- filter_conditions(csv_irf, conditions_to_plot)
csv_stat_filtered <- filter_conditions(csv_stat, conditions_to_plot)

# Calculate frequencies for each filtered dataframe
freq_homeobox <- calculate_frequency(csv_homeobox_filtered)
freq_irf <- calculate_frequency(csv_irf_filtered)
freq_stat <- calculate_frequency(csv_stat_filtered)

# Create plots for MEF and LMC conditions separately for each family

# Homeobox Family
plot_homeobox_mef <- create_plot(freq_homeobox %>% filter(Condition %in% c("I_MEF_up", "I_MEF_down", "II_MEF_up", "II_MEF_down")), "Homeobox Family Motifs [MEF]")
plot_homeobox_lmc <- create_plot(freq_homeobox %>% filter(Condition %in% c("I_LMC_up", "I_LMC_down", "II_LMC_up", "II_LMC_down")), "Homeobox Family Motifs [LMC]")

# IRF Family
plot_irf_mef <- create_plot(freq_irf %>% filter(Condition %in% c("I_MEF_up", "I_MEF_down", "II_MEF_up", "II_MEF_down")), "IRF Family Motifs [MEF]")
plot_irf_lmc <- create_plot(freq_irf %>% filter(Condition %in% c("I_LMC_up", "I_LMC_down", "II_LMC_up", "II_LMC_down")), "IRF Family Motifs [LMC]")

# STAT Family
plot_stat_mef <- create_plot(freq_stat %>% filter(Condition %in% c("I_MEF_up", "I_MEF_down", "II_MEF_up", "II_MEF_down")), "STAT Family Motifs [MEF]")
plot_stat_lmc <- create_plot(freq_stat %>% filter(Condition %in% c("I_LMC_up", "I_LMC_down", "II_LMC_up", "II_LMC_down")), "STAT Family Motifs [LMC]")

# Combine MEF and LMC plots side by side for each family

# Homeobox Family Combined Plot
combined_homeobox_plot <- plot_grid(plot_homeobox_mef, plot_homeobox_lmc, labels = c("MEF", "LMC"), ncol = 2, align = "v")

# IRF Family Combined Plot
combined_irf_plot <- plot_grid(plot_irf_mef, plot_irf_lmc, labels = c("MEF", "LMC"), ncol = 2, align = "v")

# STAT Family Combined Plot
combined_stat_plot <- plot_grid(plot_stat_mef, plot_stat_lmc, labels = c("MEF", "LMC"), ncol = 2, align = "v")

# Combine all family plots into a single figure with main title
combined_plot <- plot_grid(combined_homeobox_plot, combined_irf_plot, combined_stat_plot, labels = c("Homeobox", "IRF", "STAT"), ncol = 1, align = "v")

# Add the main title
main_title <- ggdraw() + draw_label("Analysis I and II for Homeobox, IRF, and STAT Families (MEF and LMC)", fontface = 'bold', size = 20)

# Combine the main title with the combined plots
final_plot <- plot_grid(main_title, combined_plot, ncol = 1, rel_heights = c(0.1, 1))

# Save the final plot
ggsave("Homeobox_IRF_STAT_Analysis_I_II_MEF_LMC.png", plot = final_plot, width = 18, height = 24)
