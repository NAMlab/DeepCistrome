# R Script to Load, Inspect, and Visualize the Merged Biological Data

# Install and load necessary packages
if (!require("dplyr")) install.packages("dplyr")
if (!require("stringr")) install.packages("stringr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("viridis")) install.packages("viridis") # For colorblind-friendly palettes
if (!require("ggpubr")) install.packages("ggpubr")   # For adding stats to plots
if (!require("rstatix")) install.packages("rstatix") # For easier statistical calculations
if (!require("scales")) install.packages("scales") # For color palettes
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(rstatix)
library(scales)

# --- 1. DEFINE FILE PATH ---
file_to_load <- "Peleke-etal-Suppl_table_GWAS-STARRseq_stats_sz2025h.csv"


# --- 2. LOAD THE DATASET ---
# Check if the file exists before trying to load it
if (!file.exists(file_to_load)) {
  stop(paste("Error: Cannot find the file ->", file_to_load))
}
loaded_data <- read.csv(file_to_load, header = TRUE, stringsAsFactors = FALSE)


# --- 3. DATA PREPARATION FOR VISUALIZATION ---
# Extract TFBS family and determine binding change from the '_tnt' columns.
viz_data <- loaded_data %>%
  mutate(
    TFBS_Family = str_extract(SequenceHeader, "MYB|bHLH|BZR|BES1|bZIP|WRKY")
  ) %>%
  filter(!is.na(TFBS_Family)) %>%
  # ADDED FILTER: Exclude Wild Type rows where relative strength is zero from all analyses
  filter(!(rel_strength == 0 & str_detect(SequenceHeader, "_WT"))) %>%
  # Create a new column 'Binding_Change' based on the value in the corresponding TFBS_tnt column
  mutate(
    Binding_Change_Value = case_when(
      TFBS_Family == "MYB"   ~ MYB_tnt,
      TFBS_Family == "bHLH"  ~ bHLH_tnt,
      #      TFBS_Family == "BZR"   ~ BZR_tnt,
      TFBS_Family == "BES1"  ~ BES1_tnt, # FIX: Corrected typo from BES1_nt to BES1_tnt
      TFBS_Family == "bZIP"  ~ bZIP_tnt,
      TFBS_Family == "WRKY"  ~ WRKY_tnt,
      TRUE ~ 0 # Default to 0 if no match
    ),
    Binding_Change = factor(case_when(
      Binding_Change_Value == -1  ~ "Gain",
      Binding_Change_Value == 1 ~ "Loss",
      TRUE                      ~ "No Change"
    ), levels = c("Gain", "Loss", "No Change")), # Set factor levels for legend order
    is_significant_bootstrap = as.factor(is_significant_bootstrap)
  )

# Create a separate dataframe for only the bootstrap-significant results
significant_data <- viz_data %>%
  filter(toupper(as.character(is_significant_bootstrap)) == "TRUE" & rel_strength != 0)


# --- 4. VISUALIZATION: INDIVIDUAL TFBS FAMILY VS. ALL ---
# This reusable function performs the comparison and generates a boxplot.
create_family_vs_all_plot <- function(data, title_suffix) {
  
  # --- 1. Data Preparation ---
  data_for_comparison <- data %>%
    filter(Binding_Change %in% c("Gain", "Loss"))
  
  if(nrow(data_for_comparison) == 0) {
    print(paste("Skipping 'Family vs. All' plot for", title_suffix, " - no Gain/Loss data found."))
    return(NULL)
  }
  
  all_events_data <- data_for_comparison %>%
    mutate(TFBS_Family = "All")
  
  family_levels <- c("All", sort(unique(as.character(data_for_comparison$TFBS_Family))))
  
  combined_plot_data <- bind_rows(data_for_comparison, all_events_data) %>%
    mutate(TFBS_Family = factor(TFBS_Family, levels = family_levels))
  
  # --- 2. Robust Statistical Testing (REVISED & SIMPLIFIED) ---
  # This function identifies valid groups and then lets rstatix perform all pairwise comparisons.
  perform_robust_tests <- function(data, ref_group = "All") {
    # Identify families with at least 3 data points to test
    valid_families <- data %>%
      filter(TFBS_Family != ref_group) %>%
      group_by(TFBS_Family) %>%
      summarise(n = n(), .groups = 'drop') %>%
      filter(n >= 3) %>%
      pull(TFBS_Family)
    
    if (length(valid_families) == 0) return(tibble())
    
    # Filter the data to include ONLY the reference group and the valid families
    test_data <- data %>%
      filter(TFBS_Family %in% c(ref_group, as.character(valid_families))) %>%
      droplevels() # Good practice to remove unused factor levels
    
    # Perform the Wilcoxon test once. rstatix will handle all pairwise comparisons against the ref_group.
    rstatix::wilcox_test(test_data, rel_strength ~ TFBS_Family, ref.group = ref_group) %>%
      add_significance(p.col = "p.adj", output.col = "p.adj.signif")
  }
  
  stat_test_gain <- perform_robust_tests(filter(combined_plot_data, Binding_Change == "Gain"))
  stat_test_loss <- perform_robust_tests(filter(combined_plot_data, Binding_Change == "Loss"))
  
  stat_test_combined <- bind_rows(
    if(nrow(stat_test_gain) > 0) mutate(stat_test_gain, Binding_Change = "Gain"),
    if(nrow(stat_test_loss) > 0) mutate(stat_test_loss, Binding_Change = "Loss")
  )
  
  # --- Print the statistical results table ---
  print(paste("--- Statistical Test Results for:", title_suffix, "---"))
  if (nrow(stat_test_combined) > 0) {
    print(stat_test_combined)
  } else {
    print("No valid comparisons with n>=3 found.")
  }
  
  # --- 3. Plotting ---
  # Define a consistent color palette
  family_names <- setdiff(family_levels, "All")
  family_colors <- scales::viridis_pal(option = "turbo")(length(family_names))
  names(family_colors) <- family_names
  color_palette <- c("All" = "grey70", family_colors)
  
  plot <- ggplot(combined_plot_data, aes(x = TFBS_Family, y = rel_strength, fill = TFBS_Family)) +
    geom_boxplot(show.legend = FALSE) +
    facet_wrap(~ Binding_Change, scales = "free_y") +
    scale_fill_manual(values = color_palette, drop = FALSE) +
    coord_cartesian(ylim = c(-2.5, 2.5)) +
    labs(
      title = "Effect of Individual TFBS Families vs. All Combined Events",
      subtitle = paste("Data:", title_suffix, "(Significance vs. 'All' group, Wilcoxon test, n>=3)"),
      x = "TFBS Family (vs. All)",
      y = "log2(CRE strength rel. to Col0)"
    ) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (nrow(stat_test_combined) > 0) {
    plot <- plot + stat_pvalue_manual(
      stat_test_combined, 
      label = "p.adj.signif", # FIX: Use the correct significance column name
      y.position = 2.3,
      xmin = "group2", xmax = "group2", # Use group2 which is the family name
      inherit.aes = FALSE,
      size = 6
    )
  }
  
  print(paste("--- Displaying Plot: Individual TFBS Family vs. All for", title_suffix, "---"))
  print(plot)
  
  # --- 4. Save the plot ---
  plot_filename_base <- paste0("TFBS_family_vs_all_comparison_", gsub(" ", "_", title_suffix))
  ggsave(paste0(plot_filename_base, ".png"), plot, width = 10, height = 7, units = "in", dpi = 300)
  ggsave(paste0(plot_filename_base, ".pdf"), plot, width = 10, height = 7, units = "in")
  
  return(plot)
}

# --- 5. Generate and Save the Plots ---
# Generate the plot for the full dataset
plot_full <- create_family_vs_all_plot(data = viz_data, title_suffix = "Full Dataset")

# Generate the plot for the significant data only
plot_significant <- create_family_vs_all_plot(data = significant_data, title_suffix = "Significant Bootstrap Only")

print("--- Script finished ---")

