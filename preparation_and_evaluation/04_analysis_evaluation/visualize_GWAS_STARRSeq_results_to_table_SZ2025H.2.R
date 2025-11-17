# R Script to Load, Inspect, and Visualize the Merged Biological Data

# --- 0. SETUP: Install and load necessary packages ---
# This block checks if packages are installed and loads them.
if (!require("dplyr")) install.packages("dplyr")
if (!require("stringr")) install.packages("stringr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("viridis")) install.packages("viridis") # For colorblind-friendly palettes
if (!require("ggpubr")) install.packages("ggpubr") # For adding stats to plots
if (!require("pheatmap")) install.packages("pheatmap") # For the co-occurrence heatmap
if (!require("ggrepel")) install.packages("ggrepel")
if (!require("tidyr")) install.packages("tidyr") # For pivot_longer()
if (!require("RColorBrewer")) install.packages("RColorBrewer")
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(pheatmap)
library(ggrepel)
library(tidyr)
library(RColorBrewer)
library(knitr)

# --- 1. DEFINE FILE PATH ---
# Please update this path to the location of your CSV file
file_to_load <- "~/Desktop/deepCIS_calc/studies/deepCIS/candidate_screening/Peleke-etal-Suppl_table_GWAS-STARRseq_stats_sz2025h.csv"


# --- 2. LOAD THE DATASET ---
# Check if the file exists before trying to load it
if (!file.exists(file_to_load)) {
  stop(paste("Error: Cannot find the file ->", file_to_load))
}
loaded_data <- read.csv(file_to_load, header = TRUE, stringsAsFactors = FALSE)

# --- DIAGNOSTIC STEP 1: Inspect the loaded data ---
print("--- DIAGNOSTIC: Initial Data Load ---")
print(paste("Dimensions of loaded data:", paste(dim(loaded_data), collapse = " x ")))


# --- 3. DATA PREPARATION FOR VISUALIZATION ---
# Extract TFBS family and determine binding change from the '_tnt' columns.
viz_data <- loaded_data %>%
  mutate(
    TFBS_Family = str_extract(SequenceHeader, "MYB|bHLH|BZR|BES1|bZIP|WRKY")
  ) %>%
  filter(!is.na(TFBS_Family)) %>%
  # Create a new column 'Binding_Change' based on the value in the corresponding TFBS_tnt column
  mutate(
    construct_type = str_extract(SequenceHeader, "[^_]+$"),
    Binding_Change_Value = case_when(
      TFBS_Family == "MYB"    ~ MYB_tnt,
      TFBS_Family == "bHLH"   ~ bHLH_tnt,
      TFBS_Family == "BZR"    ~ BZR_tnt,
      TFBS_Family == "BES1"   ~ BES1_tnt,
      TFBS_Family == "bZIP"   ~ bZIP_tnt,
      TFBS_Family == "WRKY"   ~ WRKY_tnt,
      TRUE ~ 0 # Default to 0 if no match
    ),
    Binding_Change = factor(case_when(
      Binding_Change_Value == -1 ~ "Gain",
      Binding_Change_Value == 1  ~ "Loss",
      TRUE                     ~ "No Change"
    ), levels = c("Gain", "Loss", "No Change")), # Set factor levels for legend order
    is_significant_bootstrap = as.factor(is_significant_bootstrap)
  )

# --- DIAGNOSTIC STEP 2: Inspect the prepared data ---
print("--- DIAGNOSTIC: Data after TFBS Extraction and Binding Change Calculation ---")
print(paste("Dimensions of viz_data:", paste(dim(viz_data), collapse = " x ")))
if (nrow(viz_data) == 0) {
  print("WARNING: viz_data is empty. TFBS_Family extraction might have failed for all rows.")
} else {
  print("Head of viz_data (showing new columns):")
  print(head(select(viz_data, SequenceHeader, construct_type, TFBS_Family, Binding_Change, rel_strength, is_significant_bootstrap)))
}

################################################################################
# --- 4. DEFINE CONSISTENT VISUAL STYLES ---
################################################################################

# Create a consistent color palette for Control vs. Variant
tf_color_palette_2 <- c("Control" = "grey59", "Variant" = "white")

# Create a consistent color palette for TFBS Families
unique_TFBS_Families <- unique(viz_data$TFBS_Family)
tf_color_palette <- setNames(
  RColorBrewer::brewer.pal(length(unique_TFBS_Families), "Set2"),
  unique_TFBS_Families
)

print("--- Consistent Color Palette Created ---")
print(tf_color_palette)


# Create a reusable custom theme for consistent, large, and readable plots
theme_custom <- theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.2)),
    plot.subtitle = element_text(hjust = 0.5, size = rel(1.1)),
    axis.text.x = element_text(angle = 60, hjust = 1, size = rel(1.2)),
    legend.position = "bottom",
    panel.grid.major = element_line(colour = "grey75", linewidth = 0.5),
    panel.grid.minor = element_line(colour = "grey75", linewidth = 0.25),
    axis.line = element_line(colour = "black")
  )

################################################################################
# --- 5. SEPARATE VISUALIZATIONS: GAIN & LOSS PLOTS ---
################################################################################

# --- Prepare a single, unified dataset for boxplots and spaghetti lines ---
# This modern approach uses pivot_longer to create a tidy dataset from the
# individual replicate data. This is more statistically sound for a paired analysis.
plot_data_long <- viz_data %>%
  # Filter for the relevant events and regions first
  filter(Binding_Change %in% c("Gain", "Loss"), cis_region %in% c("promoter", "terminator")) %>%
  # Add a unique ID for each original row to track pairs for lines and stats
  mutate(id = row_number()) %>%
  # Pivot the replicate columns (ctrl_1, var_1, ctrl_2, var_2) into a long format
  pivot_longer(
    cols = c(ctrl_1, var_1, ctrl_2, var_2),
    names_to = c("Condition", "Replicate"), # Creates 'Condition' and 'Replicate' columns
    names_sep = "_",                      # Splits "ctrl_1" into "ctrl" and "1"
    values_to = "Value"                   # Numeric values go into the 'Value' column
  ) %>%
  # Clean up the new 'Condition' column to be a factor for plotting
  mutate(
    Condition = factor(Condition,
                       levels = c("ctrl", "var"),
                       labels = c("Control", "Variant"))
  ) %>%
  # Calculate the direction of change for each pair (for coloring spaghetti lines)
  group_by(id, Replicate) %>%
  mutate(
    Direction = case_when(
      last(Value) > first(Value) ~ "Increasing",
      last(Value) < first(Value) ~ "Decreasing",
      TRUE ~ "No Change"
    )
  ) %>%
  ungroup() # Ungroup after the operation

# --- PLOT 1: GAIN EVENTS ONLY ---
# 1. Filter data and calculate sample sizes for GAIN events
plot_data_gain <- plot_data_long %>% filter(Binding_Change == "Gain")
sample_sizes_gain <- plot_data_gain %>%
  group_by(TFBS_Family) %>%
  summarise(n = n_distinct(id), .groups = 'drop')

# 2. Generate the GAIN plot
if (nrow(plot_data_gain) > 0) {
  gain_plot <- ggplot(
    plot_data_gain,
    aes(x = Condition, y = Value, fill = Condition)
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    # Draw spaghetti lines for ALL pairs, not just significant ones
    geom_line(
      aes(group = interaction(id, Replicate), color = Direction),
      alpha = 0.5,
      linewidth = 0.5
    ) +
    geom_jitter(width = 0, size = 0.7, color = "gray20") +
    # Facet by TFBS Family, arranged in one horizontal row
    facet_wrap(~ TFBS_Family, nrow = 1) +
    stat_compare_means(
      paired = TRUE, method = "wilcox.test", label = "p.format",
      label.y.npc = 0.95
    ) +
    geom_text(
      data = sample_sizes_gain,
      aes(label = paste("n =", n), x = 1.5, y = -Inf),
      vjust = -0.5, inherit.aes = FALSE, size = 4
    ) +
    scale_fill_manual(values = tf_color_palette_2) +
    scale_color_manual(
      name = "Expression Change",
      values = c("Increasing" = "red", "Decreasing" = "blue", "No Change" = "grey50")
    ) +
    labs(
      title = "Control vs. Variant Expression (TFBS Gain Events)",
      x = "Condition", y = "Delta to Baseline Expression", fill = "Condition"
    ) +
    theme_custom +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      strip.background = element_rect(fill = "grey85"),
      strip.text = element_text(face = "bold")
    )
  
  print("--- Displaying Gain Events Plot ---")
  print(gain_plot)
  
  # Save the plot
  ggsave(
    filename = "GWAS_TF_Gain_Events_Horizontal.pdf", plot = gain_plot,
    width = 12, height = 6, units = "in", dpi = 600
  )
  ggsave(
    filename = "GWAS_TF_Gain_Events_Horizontal.svg", plot = gain_plot,
    width = 12, height = 6, units = "in", dpi = 600
  )
} else {
  print("Skipping Gain plot: No data available.")
}

# --- PLOT 2: LOSS EVENTS ONLY ---
# 1. Filter data and calculate sample sizes for LOSS events
plot_data_loss <- plot_data_long %>% filter(Binding_Change == "Loss")
sample_sizes_loss <- plot_data_loss %>%
  group_by(TFBS_Family) %>%
  summarise(n = n_distinct(id), .groups = 'drop')

# 2. Generate the LOSS plot
if (nrow(plot_data_loss) > 0) {
  loss_plot <- ggplot(
    plot_data_loss,
    aes(x = Condition, y = Value, fill = Condition)
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    # Draw spaghetti lines for ALL pairs
    geom_line(
      aes(group = interaction(id, Replicate), color = Direction),
      alpha = 0.5,
      linewidth = 0.5
    ) +
    geom_jitter(width = 0, size = 0.7, color = "gray20") +
    # Facet by TFBS Family, arranged in one horizontal row
    facet_wrap(~ TFBS_Family, nrow = 1) +
    stat_compare_means(
      paired = TRUE, method = "wilcox.test", label = "p.format",
      label.y.npc = 0.95
    ) +
    geom_text(
      data = sample_sizes_loss,
      aes(label = paste("n =", n), x = 1.5, y = -Inf),
      vjust = -0.5, inherit.aes = FALSE, size = 4
    ) +
    scale_fill_manual(values = tf_color_palette_2) +
    scale_color_manual(
      name = "Expression Change",
      values = c("Increasing" = "red", "Decreasing" = "blue", "No Change" = "grey50")
    ) +
    labs(
      title = "Control vs. Variant Expression (TFBS Loss Events)",
      x = "Condition", y = "Delta to Baseline Expression", fill = "Condition"
    ) +
    theme_custom +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      strip.background = element_rect(fill = "grey85"),
      strip.text = element_text(face = "bold")
    )
  
  print("--- Displaying Loss Events Plot ---")
  print(loss_plot)
  
  # Save the plot
  ggsave(
    filename = "GWAS_TF_Loss_Events_Horizontal.pdf", plot = loss_plot,
    width = 12, height = 6, units = "in", dpi = 600
  )
  ggsave(
    filename = "GWAS_TF_Loss_Events_Horizontal.svg", plot = loss_plot,
    width = 12, height = 6, units = "in", dpi = 600
  )
} else {
  print("Skipping Loss plot: No data available.")
}


################################################################################
# --- DATA DEDUPLICATION & SUBSEQUENT PLOTS ---
################################################################################

# Remove identical SNPs (due to overlapping prediction windows) by keeping the one with the max effect
viz_data_deduplicated <- viz_data %>%
  group_by(SNP_position) %>%
  slice_max(order_by = abs(rel_strength), n = 1, with_ties = FALSE) %>%
  ungroup()

# --- PLOT 3: DEEP DIVE BAR PLOT FOR A SINGLE PHENOTYPE ---
target_phenotype_deep_dive <- "DTFplantingSummerLocSweden2009"

phenotype_specific_data_unfiltered <- viz_data_deduplicated %>%
  filter(
    study.phenotype.name == target_phenotype_deep_dive,
    str_detect(SequenceHeader, "_mutA") # This keeps only the mutant variants
  ) %>%
  mutate(xaxis_label = paste(SNP_position, gene_id, Binding_Change, sep = "\n"))

if(nrow(phenotype_specific_data_unfiltered) > 0) {
  plot3_deep_dive <- ggplot(phenotype_specific_data_unfiltered,
                            aes(x = reorder(xaxis_label, rel_strength),
                                y = rel_strength,
                                fill = TFBS_Family)) +
    geom_col(width = 0.6, alpha = 0.8) +
    geom_hline(yintercept = 0, color = "grey50") +
    scale_fill_manual(values = tf_color_palette) +
    labs(
      title = paste("Days of Flowering Phenotype"),
      x = "SNP Position, Gene ID, and Binding Status",
      y = "log2(CRE strength rel. to Col0)",
      fill = "TFBS Family"
    ) +
    theme_custom +
    scale_y_continuous(breaks = seq(-2, 2, by = 0.25)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14))
  
  print("--- Displaying Single Phenotype Deep Dive Bar Plot ---")
  print(plot3_deep_dive)
  
  # Save the plot now that it has been created
  ggsave(
    filename = "GWAS_phenotypes_DTF.pdf", plot = plot3_deep_dive,
    width = 6, height = 6, units = "in", dpi = 600
  )
  ggsave(
    filename = "GWAS_phenotypes_DTF.svg", plot = plot3_deep_dive,
    width = 6, height = 6, units = "in", dpi = 600
  )
  
} else {
  print(paste("Skipping plot: no data found for phenotype:", target_phenotype_deep_dive))
}


# --- PLOT 4: COMBINED BAR PLOT FOR BACTERIAL RESISTANCE PHENOTYPES ---
combined_avr_data <- viz_data_deduplicated %>%
  filter(
    # *** CORRECTION: Using str_detect with '^' anchor for better compatibility ***
    str_detect(study.phenotype.name, "^avr"),
    str_detect(SequenceHeader, "_mutA") # Keeps only mutant variants
  ) %>%
  mutate(xaxis_label = paste(study.phenotype.name, SNP_position, gene_id, Binding_Change, sep = "\n"))


if (nrow(combined_avr_data) > 0) {
  print(paste("Found", nrow(combined_avr_data), "data points across all 'avr' phenotypes. Generating combined plot."))
  
  combined_avr_plot <- ggplot(combined_avr_data,
                              aes(x = reorder(xaxis_label, rel_strength),
                                  y = rel_strength,
                                  fill = TFBS_Family)) +
    geom_col(width = 0.7, alpha = 0.8) +
    geom_hline(yintercept = 0, color = "grey50") +
    scale_fill_manual(values = tf_color_palette) +
    labs(
      title = "Bacterial Disease Resistance Phenotypes",
      x = "Phenotype, SNP Position, Gene ID, and Binding Status",
      y = "log2(CRE strength rel. to Col0)",
      fill = "TFBS Family"
    ) +
    theme_custom +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14))
  
  print(combined_avr_plot)
  
  # Save the plot now that it has been created
  ggsave(
    filename = "GWAS_phenotypes_RESIST.pdf", plot = combined_avr_plot,
    width = 6, height = 6, units = "in", dpi = 600
  )
  ggsave(
    filename = "GWAS_phenotypes_RESIST.svg", plot = combined_avr_plot,
    width = 6, height = 6, units = "in", dpi = 600
  )
  
} else {
  print("No data found for any phenotypes starting with 'avr'. No plot was generated.")
}

################################################################################
# --- 7. GENERATE SUMMARY STATISTICS & TEST TABLES ---
################################################################################
# --- Statistical Test Results Table ---
# This performs the paired Wilcoxon test for each TFBS Family and Binding_Change
# combination, exactly matching the tests shown on the plots.
stat_tests <- plot_data_long %>%
  # Group by the facets used in the plots
  group_by(TFBS_Family, Binding_Change) %>%
  # group_modify is a clean way to apply a function to each group subset
  group_modify(~ {
    # --- START: NEW & CORRECTED LOGIC ---
    # Reshape the data for this group to create explicit 'Control' and 'Variant' columns
    wide_data <- .x %>%
      pivot_wider(
        id_cols = c(id, Replicate),
        names_from = Condition,
        values_from = Value
      )
    
    # Remove any rows that are not complete pairs (e.g., a Control without a Variant)
    complete_pairs <- na.omit(wide_data)
    
    # Wilcoxon test needs at least a few pairs to run meaningfully.
    # A minimum of ~6 pairs is a reasonable threshold.
    if (nrow(complete_pairs) < 6) {
      return(tibble(statistic = NA, p.value = NA))
    }
    
    # Perform the Wilcoxon test on the correctly paired data
    test_result <- wilcox.test(
      complete_pairs$Control,
      complete_pairs$Variant,
      paired = TRUE,
      exact = FALSE # Set to FALSE for larger samples
    )
    # --- END: NEW & CORRECTED LOGIC ---
    
    # Return a tidy data frame of the results using broom::tidy
    broom::tidy(test_result)
  }) %>%
  ungroup() %>%
  # This select call will now work correctly for all rows
  select(
    TFBS_Family,
    Binding_Change,
    v_statistic = statistic,
    p_value = p.value
  )

print(stat_tests)
