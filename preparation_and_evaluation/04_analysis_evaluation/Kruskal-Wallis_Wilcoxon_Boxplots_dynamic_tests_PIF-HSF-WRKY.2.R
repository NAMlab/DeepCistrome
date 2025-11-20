# --- 0. Setup: Load Libraries ---
library(dplyr)
library(stringr)
library(ggplot2)
library(broom)
library(ggpubr)
library(tidyr) # Added for pivot_longer

cat("Current working directory:", getwd(), "\n")

# --- 1. Load and Perform Initial Data Cleaning ---
file_path <- "TFBS_mutated_dCIS_candidate_predictions_condition.csv"
if (!file.exists(file_path)) {
  stop(paste("Error: File not found at", file_path))
}
df <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
cat("Data loaded successfully. Dimensions:", dim(df)[1], "rows,", dim(df)[2], "columns.\n")

# --- Essential Column Checks and Processing ---
if (!("condition" %in% colnames(df))) {
  stop("Error: The crucial 'condition' column was not found in the data file.")
}
first_tfbs_pred_col_name <- "ABI3VP1"
if (!(first_tfbs_pred_col_name %in% colnames(df))) {
  stop(paste("Error: Expected first TFBS predictor column '", first_tfbs_pred_col_name, "' not found.", sep=""))
}
tfbs_pred_cols_start_idx <- which(colnames(df) == first_tfbs_pred_col_name)
original_predictor_column_names <- colnames(df)[tfbs_pred_cols_start_idx:ncol(df)]

# Convert TFBS predictor columns to integers
for (col in original_predictor_column_names) {
  if (col %in% colnames(df)) {
    if (is.character(df[[col]])) { df[[col]] <- gsub("00", "0", df[[col]]) }
    df[[col]] <- suppressWarnings(as.integer(df[[col]]))
  }
}
if (!("rel_strength" %in% colnames(df))) stop("Error: 'rel_strength' column not found.")
df$rel_strength <- as.numeric(as.character(df$rel_strength))
cat("'rel_strength' and TFBS predictor columns processed.\n")


# --- 2. Create Core Identifier Columns ---
df <- df %>%
  mutate(
    Species = substr(SequenceHeader, 1, 2),
    TargetedTFBS_from_header = sapply(SequenceHeader, function(h) if(length(p <- unlist(str_split(h, "_"))) >= 2) p[length(p)-1] else NA_character_),
    gene_id = sapply(strsplit(as.character(SequenceHeader), "_"), function(x) if(length(x) >= 2) x[2] else NA_character_),
    SequenceType = ifelse(endsWith(as.character(SequenceHeader), "_WT"), "WT", "Mutant")
  ) %>%
  filter(!is.na(gene_id), !is.na(TargetedTFBS_from_header), !is.na(rel_strength))

if(nrow(df) == 0) stop("No data after initial parsing of SequenceHeader.")

df <- df %>% mutate(ConstructID = paste(gene_id, TargetedTFBS_from_header, sep = "_"))
cat("Core identifier columns created.\n")


# --- 3. Define the Data Processing Function ---
# This function now processes data for a given analysis and returns a tidy dataframe,
# but it does NOT plot. Plotting is handled later.
process_tf_data <- function(base_df, target_tf) {
  
  target_df_initial <- base_df %>% filter(TargetedTFBS_from_header == target_tf)
  if(nrow(target_df_initial) == 0) {
    return(NULL) # Return NULL if no data
  }
  
  wt_rel_strength <- target_df_initial %>%
    filter(SequenceType == "WT") %>%
    group_by(ConstructID) %>%
    summarise(WT_rel_strength = mean(rel_strength, na.rm = TRUE), .groups = 'drop') %>%
    filter(!is.na(WT_rel_strength))
  
  if(nrow(wt_rel_strength) == 0) {
    return(NULL)
  }
  
  df_with_deltas <- target_df_initial %>%
    left_join(wt_rel_strength, by = "ConstructID") %>%
    mutate(delta_rel_strength = ifelse(SequenceType == "Mutant", rel_strength - WT_rel_strength, 0)) %>%
    filter(!is.na(WT_rel_strength))
  
  tf_cols_for_binding_change_all <- intersect(original_predictor_column_names, colnames(df_with_deltas))
  df_with_changes <- df_with_deltas
  
  if(length(tf_cols_for_binding_change_all) > 0) {
    safe_first <- function(x) {
      val <- na.omit(x)
      if (length(val) == 0) return(NA_integer_)
      return(first(val))
    }
    
    wt_binding_status_all_tfs <- df_with_deltas %>%
      filter(SequenceType == "WT") %>%
      group_by(ConstructID) %>%
      summarise(across(all_of(tf_cols_for_binding_change_all), safe_first), .groups = 'drop')
    
    if(nrow(wt_binding_status_all_tfs) > 0){
      colnames(wt_binding_status_all_tfs)[-1] <- paste0("WT_", colnames(wt_binding_status_all_tfs)[-1])
      df_with_changes <- df_with_deltas %>% left_join(wt_binding_status_all_tfs, by = "ConstructID")
      
      for (tf_col in tf_cols_for_binding_change_all) {
        wt_status_col_name <- paste0("WT_", tf_col)
        binding_change_col_name <- paste0(tf_col, "_Change")
        if (tf_col %in% colnames(df_with_changes) && wt_status_col_name %in% colnames(df_with_changes)) {
          df_with_changes <- df_with_changes %>%
            mutate(!!binding_change_col_name := case_when(
              SequenceType == "Mutant" & .data[[wt_status_col_name]] == 1 & .data[[tf_col]] == 0 ~ "Loss",
              SequenceType == "Mutant" & .data[[wt_status_col_name]] == 0 & .data[[tf_col]] == 1 ~ "Gain",
              SequenceType == "Mutant" & .data[[wt_status_col_name]] == 1 & .data[[tf_col]] == 1 ~ "Retained",
              SequenceType == "Mutant" & .data[[wt_status_col_name]] == 0 & .data[[tf_col]] == 0 ~ "Unbound",
              TRUE ~ NA_character_
            ))
        }
      }
    }
  }
  
  df_stat_input <- df_with_changes %>% filter(SequenceType == "Mutant", !is.na(delta_rel_strength))
  if (nrow(df_stat_input) < 10) { 
    return(NULL)
  }
  
  plot_data_long <- df_stat_input %>%
    select(condition, ConstructID, rel_strength, ends_with("_Change")) %>%
    pivot_longer(
      cols = ends_with("_Change"),
      names_to = "TFBS_Family",
      values_to = "Scenario_Raw",
      names_pattern = "(.*)_Change"
    ) %>%
    filter(!is.na(Scenario_Raw)) %>%
    mutate(Scenario = recode_factor(as.character(Scenario_Raw),
                                    "Loss"     = "Predicted Loss",
                                    "Gain"     = "Predicted Gain",
                                    "Retained" = "Predicted Retained",
                                    "Unbound"  = "Consistently Unbound",
                                    .default = NA_character_
    )) %>%
    filter(!is.na(Scenario))
  
  return(plot_data_long)
}

# --- 4. Run Analyses, Filter by Significance, and Plot ---

# Define the analyses you want to run
analyses_to_run <- list(
  list(target = "PIF", partner = "bHLH"),
  list(target = "HSF", partner = "HSF"),
  list(target = "WRKY", partner = "WRKY")
)

# Loop through each analysis type
for (analysis in analyses_to_run) {
  target_tf_name <- analysis$target
  partner_tf_name <- analysis$partner
  
  cat(paste0("\n\n############################################################\n"))
  cat(paste0("###   PROCESSING AND PLOTTING FOR TARGET: ", toupper(target_tf_name), "   ###\n"))
  cat(paste0("############################################################\n"))
  
  # 1. Collect processed data from all conditions for the current analysis
  all_conditions_data <- process_tf_data(df, target_tf_name)
  
  if (is.null(all_conditions_data) || nrow(all_conditions_data) == 0) {
    cat(paste("No data found for", target_tf_name, "across all conditions. Skipping.\n"))
    next
  }
  
  # 2. Perform Statistical Filtering
  cat("--- Filtering for statistically significant factors (p < 0.05) ---\n")
  
  # Calculate Kruskal-Wallis p-values for each group
  kruskal_results <- all_conditions_data %>%
    group_by(condition, TFBS_Family) %>%
    filter(n_distinct(Scenario) > 1) %>% # Kruskal test needs more than 1 group
    summarise(kruskal_p = kruskal.test(rel_strength ~ Scenario)$p.value, .groups = "drop")
  
  # Calculate pairwise Wilcoxon p-values for each group
  # Note: This can be slow on very large datasets
  pairwise_results <- all_conditions_data %>%
    group_by(condition, TFBS_Family) %>%
    filter(n_distinct(Scenario) > 1) %>%
    # Use nest/map approach for robustness
    nest() %>%
    mutate(
      pairwise_p = map(data, ~ broom::tidy(pairwise.wilcox.test(.$rel_strength, .$Scenario, p.adjust.method = "none")))
    ) %>%
    unnest(pairwise_p) %>%
    select(condition, TFBS_Family, p.value) %>%
    rename(wilcox_p = p.value)
  
  # Identify groups with at least one significant pairwise comparison
  significant_pairwise <- pairwise_results %>%
    filter(wilcox_p < 0.05) %>%
    select(condition, TFBS_Family) %>%
    distinct()
  
  # Identify groups with a significant overall Kruskal-Wallis test
  significant_kruskal <- kruskal_results %>%
    filter(kruskal_p < 0.05) %>%
    select(condition, TFBS_Family) %>%
    distinct()
  
  # Find the intersection: groups that are significant by BOTH criteria
  final_significant_groups <- inner_join(significant_kruskal, significant_pairwise, by = c("condition", "TFBS_Family"))
  
  if(nrow(final_significant_groups) == 0) {
    cat(paste("No factors met the significance criteria for", target_tf_name, "analysis. Skipping plot.\n"))
    next
  }
  
  cat(paste("Found", nrow(final_significant_groups), "significant factor-condition combinations to plot.\n"))
  
  # Filter the original data to keep only these significant groups for plotting
  final_plot_data <- all_conditions_data %>%
    inner_join(final_significant_groups, by = c("condition", "TFBS_Family")) %>%
    # Create new facet label with asterisk for KW significance
    mutate(TFBS_Family_Label = paste(TFBS_Family, "*"))
  
  # 3. Generate the Comparative Plot
  
  # Define all pairwise comparisons to be made for the Wilcoxon test
  my_comparisons <- combn(levels(final_plot_data$Scenario), 2, simplify = FALSE)
  
  comparative_plot <- ggplot(final_plot_data, aes(x = Scenario, y = rel_strength, fill = Scenario)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, show.legend = FALSE) + # Hide legend for fill
    geom_jitter(width = 0.2, alpha = 0.4, size = 1.2, show.legend = FALSE) +
    # Facet by the new label
    facet_grid(condition ~ TFBS_Family_Label, scales = "free_y") +
    
    # --- CHANGE HERE: Hide non-significant comparisons and reduce text size ---
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", size = 3, hide.ns = FALSE) +
    
    labs(
      title = paste("Significant TFBS Scenario Comparisons for", target_tf_name, "Constructs"),
      subtitle = "Showing factors with at least one Wilcoxon p < 0.05. '*' in title indicates overall Kruskal-Wallis p < 0.05.",
      x = NULL, 
      y = "Relative Strength of Mutant Construct"
    ) +
    # --- CHANGE HERE: Use a smaller base size and adjust individual text elements ---
    theme_bw(base_size = 16) + # Reduced base size for all text
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16), # Smaller x-axis text
      axis.title.y = element_text(size = 16), # Smaller y-axis title
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9), # Smaller main title
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 7), # Smaller subtitle
      strip.text = element_text(face = "bold", size = 16), # Smaller facet titles
      strip.background = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))+
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed", linewidth = 0.4)
  
  print(comparative_plot)
  
  # Save the plot
  results_folder_name <- "dCIS_prediction_results/filtered_comparative_plots/"
  if (!dir.exists(results_folder_name)) {
    dir.create(results_folder_name, recursive = TRUE)
  }
  plot_filename <- file.path(results_folder_name, paste0("Significant_Factors_", target_tf_name, "_Analysis.png"))
  ggsave(filename = plot_filename, plot = comparative_plot, width = max(10, 2.5 * n_distinct(final_plot_data$TFBS_Family)), height = max(8, 3 * n_distinct(final_plot_data$condition)), dpi = 300)
  cat("Saved comparative plot to:", plot_filename, "\n")
}

cat("\n\n--- All Requested Analyses are Complete ---\n")
