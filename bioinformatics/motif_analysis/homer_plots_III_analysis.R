#setwd('Desktop/Bioinformatics_projects/Bulk_RNA_seqdata_sp3/cluster_data/Homer_results/Homer_motifs_known_shortlisted/')

# Load the necessary libraries
library(ggplot2)
library(gridExtra)
library(cowplot)
library(dplyr)
library(stringr)

# Define a mapping of old condition names to the new labels
condition_labels <- c(
  "III_MEF_up" = "Higher in Sp3 [MEF]",
  "III_MEF_down" = "Lower in Sp3 [MEF]",
  "III_LMC_up" = "Higher in Sp3 [LMC]",
  "III_LMC_down" = "Lower in Sp3 [LMC]"
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
      "III_MEF_up" = "#96cb01", 
      "III_MEF_down" = "#fe9388", 
      "III_LMC_up" = "#96cb01", 
      "III_LMC_down" = "#fe9388"
    ), labels = c(
      "III_MEF_up" = "higher in Sp3", 
      "III_MEF_down" = "lower in Sp3", 
      "III_LMC_up" = "higher in Sp3", 
      "III_LMC_down" = "lower in Sp3"
    )) +
    expand_limits(y = 0)  # Ensure bars start from the axis
}

# Read the CSV files
csv4 <- read.csv('all_ZF_family_motifs_combined_homeranalysis.csv')
csv5 <- read.csv('all_homeobox_family_motifs_combined_homeranalysis.csv')
csv6 <- read.csv('all_NR_family_motifs_combined_homeranalysis.csv')

# Extract conditions
csv4 <- extract_condition(csv4)
csv5 <- extract_condition(csv5)
csv6 <- extract_condition(csv6)

# Filter dataframes for specified conditions (III_MEF and III_LMC)
csv4_filtered <- filter_conditions(csv4, c("III_MEF_up", "III_MEF_down", "III_LMC_up", "III_LMC_down"))
csv5_filtered <- filter_conditions(csv5, c("III_MEF_up", "III_MEF_down", "III_LMC_up", "III_LMC_down"))
csv6_filtered <- filter_conditions(csv6, c("III_MEF_up", "III_MEF_down", "III_LMC_up", "III_LMC_down"))

# Calculate frequencies for each filtered dataframe
freq_zf <- calculate_frequency(csv4_filtered)
freq_homeobox <- calculate_frequency(csv5_filtered)
freq_nr <- calculate_frequency(csv6_filtered)

# Ensure only the specified conditions are included
freq_zf <- freq_zf %>% filter(Condition %in% c("III_MEF_up", "III_MEF_down", "III_LMC_up", "III_LMC_down"))
freq_homeobox <- freq_homeobox %>% filter(Condition %in% c("III_MEF_up", "III_MEF_down", "III_LMC_up", "III_LMC_down"))
freq_nr <- freq_nr %>% filter(Condition %in% c("III_MEF_up", "III_MEF_down", "III_LMC_up", "III_LMC_down"))

# Create plots for MEF and LMC conditions with updated labels
plot_zf_mef <- create_plot(freq_zf %>% filter(Condition %in% c("III_MEF_up", "III_MEF_down")), "ZF Family Motifs")
plot_zf_lmc <- create_plot(freq_zf %>% filter(Condition %in% c("III_LMC_up", "III_LMC_down")), "ZF Family Motifs")
plot_homeobox_mef <- create_plot(freq_homeobox %>% filter(Condition %in% c("III_MEF_up", "III_MEF_down")), "Homeobox Family Motifs")
plot_homeobox_lmc <- create_plot(freq_homeobox %>% filter(Condition %in% c("III_LMC_up", "III_LMC_down")), "Homeobox Family Motifs")
plot_nr_mef <- create_plot(freq_nr %>% filter(Condition %in% c("III_MEF_up", "III_MEF_down")), "NR Family Motifs")
plot_nr_lmc <- create_plot(freq_nr %>% filter(Condition %in% c("III_LMC_up", "III_LMC_down")), "NR Family Motifs")

# Combine ZF, Homeobox, and NR plots into a single figure with main title
combined_plot <- plot_grid(plot_zf_mef, plot_zf_lmc, plot_homeobox_mef, plot_homeobox_lmc, plot_nr_mef, plot_nr_lmc, labels = c("A", "B", "C", "D", "E", "F"), ncol = 2, align = "v")

# Add the main title
main_title <- ggdraw() + draw_label("LMC and MEF", fontface = 'bold', size = 20)

# Combine the main title with the combined plots
final_plot <- plot_grid(main_title, combined_plot, ncol = 1, rel_heights = c(0.1, 1))

# Save the final plot
ggsave("Combined_ZF_Homeobox_NR_Family_Motifs_LMC_MEF.png", plot = final_plot, width = 16, height = 18)
