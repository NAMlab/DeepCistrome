# TFBS Consolidated Co-occurrence Network Analysis
#
# This script loads a reference dataset and multiple variant datasets.
# It builds a SINGLE, consolidated co-occurrence network where:
#  - Nodes are TFBSs.
#  - An edge exists if a co-occurrence is present in ANY variant condition.
#  - Edge thickness represents the MAXIMUM strength of co-occurrence.
#  - Edge color indicates if the most significant change (p < 0.001) was a
#    gain (red) or loss (blue) of co-occurrence.
#
# The result is one comprehensive network graph created with ggraph,
# filtered to show only nodes involved in significant interactions.

# # --- 1. SETUP: Install and Load Packages ---
# 
# # Install packages if they are not already installed
# if (!requireNamespace("igraph", quietly = TRUE)) {
#   install.packages("igraph")
# }
# if (!requireNamespace("readr", quietly = TRUE)) {
#   install.packages("readr")
# }
# if (!requireNamespace("dplyr", quietly = TRUE)) {
#   install.packages("dplyr")
# }
# if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
#   install.packages("RColorBrewer")
# }
# # Add ggraph for advanced plotting
# if (!requireNamespace("ggraph", quietly = TRUE)) {
#   install.packages("ggraph")
# }
# # Add ggplot2 for saving
# if (!requireNamespace("ggplot2", quietly = TRUE)) {
#   install.packages("ggplot2")
# }
# 

# Load the required libraries
library(igraph)
library(readr)
library(dplyr)
library(RColorBrewer)
library(ggraph)
library(ggplot2)

# --- 2. DEFINE FILE PATHS AND CONSTANTS ---

# Define the path to the data files relative to the working directory.
data_path <- "./../studies/deepCIS/genome_ranges/in_silico/dCIS_TFBS_bound_candidates_AP2-WRKY-Gbox/"

# Define the common file name components
file_prefix <- "dCIS_TFBS_bound_candidates.chr1.CACGTG"
file_suffix <- ".center250w500sq.csv"

# List of the variant identifiers
variants <- c("_+1A", "_-1A", "_+1C", "_-1C", "_+1G", "_-1G", "_+1T", "_-1T", "_-1A+1T", "_-1T+1A", "_-1C+1G", "_-1G+1C") # , 

#variants <- c("_+1A","_-1A", "_+1C", "_-1C", "_+1G", "_-1G", "_+1T", "_-1T")
#("_+1A","_-1A", "_+1C", "_-1C", "_+1G", "_-1G", "_+1T", "_-1T")
#variants <- c("_anti")
#variants <- c("_global-anti")
#variants <- c("_ref_rc")


TFBS <- "AP2EREBP"
# Create the full file paths for the variant files
variant_files <- paste0(data_path, file_prefix, variants, file_suffix)
names(variant_files) <- gsub("_", "", variants)

# Define the path for the reference file
ref_file <- paste0(data_path, file_prefix, "_ref", file_suffix)


# --- 3. LOAD AND PROCESS DATA ---

# Function to read and process a single TFBS data file
load_tfbs_data <- function(filepath) {
  if (!file.exists(filepath)) {
    stop("File not found: ", filepath)
  }
  df <- read_tsv(filepath, col_types = cols(), name_repair = "minimal")
  # Store sequence count before modifying the dataframe
  sequence_count <- nrow(df)
  df <- df %>%
    tibble::column_to_rownames(var = colnames(df)[1])
  df <- as.data.frame(sapply(df, as.numeric))
  
  # Return both the dataframe and the sequence count
  return(list(data = df, count = sequence_count))
}

# Function to calculate co-occurrence matrix
calculate_cooccurrence <- function(df) {
  m <- as.matrix(df)
  m[is.na(m)] <- 0
  # Ensure only TFBSs with some binding are included
  m <- m[, colSums(m) > 0]
  # Calculate co-occurrence
  co_occurrence_matrix <- t(m) %*% m
  diag(co_occurrence_matrix) <- 0 # No self-loops
  return(co_occurrence_matrix)
}

# Load the reference data ONCE
cat("Loading reference data from:", ref_file, "\n")
ref_info <- load_tfbs_data(ref_file)
ref_data <- ref_info$data
ref_seq_count <- ref_info$count
ref_cooc <- calculate_cooccurrence(ref_data)
cat("Reference data loaded.\n\n")


# --- 4. CALCULATE AND CONSOLIDATE MATRICES ---

cooc_matrices <- list()
seq_counts <- list()
cat("Calculating co-occurrence for all variants...\n")
for (variant_name in names(variant_files)) {
  variant_filepath <- variant_files[variant_name]
  variant_info <- load_tfbs_data(variant_filepath)
  cooc_matrices[[variant_name]] <- calculate_cooccurrence(variant_info$data)
  seq_counts[[variant_name]] <- variant_info$count
}

# Get a list of all TFBSs present across all datasets
all_tfbs_names <- unique(c(rownames(ref_cooc), unlist(lapply(cooc_matrices, rownames))))

# Create empty master matrices
master_weight_matrix <- matrix(0, nrow=length(all_tfbs_names), ncol=length(all_tfbs_names), dimnames=list(all_tfbs_names, all_tfbs_names))
master_sig_type_matrix <- matrix(NA, nrow=length(all_tfbs_names), ncol=length(all_tfbs_names), dimnames=list(all_tfbs_names, all_tfbs_names))

# Loop through each pair of TFBSs to find the max co-occurrence and most significant change
for (i in 1:nrow(master_weight_matrix)) {
  for (j in i:ncol(master_weight_matrix)) {
    if (i == j) next
    
    tfbs1 <- rownames(master_weight_matrix)[i]
    tfbs2 <- colnames(master_weight_matrix)[j]
    
    max_cooc <- -1
    min_p_value <- 1
    direction_of_sig_change <- NA
    
    ref_cooc_count <- ifelse(tfbs1 %in% rownames(ref_cooc) && tfbs2 %in% colnames(ref_cooc), ref_cooc[tfbs1, tfbs2], 0)
    
    for (variant_name in names(cooc_matrices)) {
      mat <- cooc_matrices[[variant_name]]
      variant_cooc_count <- ifelse(tfbs1 %in% rownames(mat) && tfbs2 %in% colnames(mat), mat[tfbs1, tfbs2], 0)
      
      # Update max co-occurrence for edge weight
      if (variant_cooc_count > max_cooc) max_cooc <- variant_cooc_count
      
      # Perform Fisher's test for significance of change
      contingency_table <- matrix(c(
        variant_cooc_count, seq_counts[[variant_name]] - variant_cooc_count,
        ref_cooc_count,     ref_seq_count - ref_cooc_count
      ), nrow = 2, byrow = TRUE)
      
      p_val <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      
      if (p_val < min_p_value) {
        min_p_value <- p_val
        # Check direction of change
        if(variant_cooc_count > ref_cooc_count){
          direction_of_sig_change <- "Positive"
        } else if (variant_cooc_count < ref_cooc_count){
          direction_of_sig_change <- "Negative"
        } else {
          direction_of_sig_change <- "No Change"
        }
      }
    }
    
    if (max_cooc > 0) {
      master_weight_matrix[i, j] <- master_weight_matrix[j, i] <- max_cooc
      if (min_p_value < 0.001 && direction_of_sig_change != "No Change") {
        master_sig_type_matrix[i, j] <- master_sig_type_matrix[j, i] <- direction_of_sig_change
      }
    }
  }
}

# --- 5. GENERATE THE CONSOLIDATED NETWORK GRAPH WITH GGRAPH ---

cat("\nGenerating the consolidated network graph with ggraph...\n")

# Create the full graph first
full_graph <- graph_from_adjacency_matrix(master_weight_matrix, mode = "upper", weighted = TRUE, diag = FALSE)

# Identify nodes involved in significant edges
sig_edges <- which(!is.na(master_sig_type_matrix))
sig_nodes_indices <- unique(as.vector(arrayInd(sig_edges, dim(master_sig_type_matrix))))
sig_nodes_names <- rownames(master_sig_type_matrix)[sig_nodes_indices]


if(length(sig_nodes_names) == 0){
  stop("No significant edges (p < 0.001) were found. Cannot create filtered network.")
}

# Create a subgraph containing only the significant nodes and their edges
graph <- induced_subgraph(full_graph, vids = sig_nodes_names)
graph <- delete.vertices(graph, which(degree(graph) < 1)) # Clean up any isolated nodes after filtering

if (gsize(graph) == 0) {
  stop("Analysis finished, but no connected significant edges were found to build a network.")
}

# Re-get the significant change type for the filtered graph
edge_list <- as_edgelist(graph, names = TRUE)
edge_change_type <- sapply(1:nrow(edge_list), function(i) master_sig_type_matrix[edge_list[i,1], edge_list[i,2]])
edge_change_type[is.na(edge_change_type)] <- "No Sig. Change" 

# Set edge attribute for coloring
E(graph)$change_type <- edge_change_type

# --- Setup Colors and Node Attributes ---
# MODIFIED: Explicitly create named vectors for colors to ensure correct mapping.
final_edge_colors <- c("Positive" = "orchid3", "Negative" = "cyan4", "No Sig. Change" = "grey60")

# Set up node colors. The main TFBS will be darkorange, others will be grey.
node_color_values <- c("darkorange", "grey40")
names(node_color_values) <- c(TFBS, "Other TFBS")
V(graph)$node_type <- ifelse(V(graph)$name == TFBS, TFBS, "Other TFBS")

# --- Create the Plot ---
network_plot <- ggraph(graph, layout = "fr") + 
  geom_edge_link(aes(width = weight, color = change_type), alpha = 0.8) +
  geom_node_point(aes(color = node_type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, size = 7, family = "sans") +
  # Use the explicitly created named vectors for colors
  scale_color_manual(values = node_color_values, name = "Node Type") +
  scale_edge_color_manual(values = final_edge_colors, name = "Direction of Sig. Change (p < 0.001)") +
  scale_edge_width_continuous(name = "Max Co-occurrence", range = c(0.5, 4)) +
  theme_void() +
  guides(color = guide_legend(override.aes = list(size=5))) + # Make node legend points larger
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  ggtitle("Consolidated TFBS Co-occurrence Network (Filtered for Significance)")

# --- 6. SAVE THE PLOT ---
# Ensure variant names are combined safely for the filename
safe_variant_string <- paste(gsub("[^A-Za-z0-9]", "", variants), collapse="-")
output_filename <- paste0(data_path, file_prefix,"_", safe_variant_string, ".consolidated_network.png")

cat(paste("  -> Saving consolidated network to ->", output_filename, "\n"))

ggsave(filename = output_filename, plot = network_plot, width = 14, height = 12, dpi = 300)

cat("\nNetwork analysis complete.\n")
