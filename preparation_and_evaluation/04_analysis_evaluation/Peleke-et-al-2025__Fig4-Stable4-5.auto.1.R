#This tool test the enrichment of EPMs to annotated genome features#
#SMZ 2024-09-01
# Load the necessary packages
library(dplyr)
library(tidyverse)
library(tidyr)
library(stringr)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(ggrepel)
library(igraph)
library(pheatmap)
library(ggraph)
library(ggplot2)
library(RColorBrewer)

########################################
output_file_path <- "IPM_TFBS_summary_statistics.csv"
# Set the working directory to the location of the BED file
setwd("/home/ibg-4/Desktop/deepCIS_calc/")
# Define the path to your CSV file
file_path <- "./studies/deepCIS/ipmArthDAPmotifs/extracted_sites_predictions_chrom_1.csv"
# Define the file that contains similar ids
file_path11 <- "./studies/deepCIS/ipmArthDAPmotifs/extracted_ranges_1.txt"
################################################################
file_path22 <- "./studies/deepCIS/ipmArthDAPmotifs/predicted_percentages_per_family.csv"
# Import the CSV file with accuracies
predicted_percentages <- read.csv(file_path22)
head(predicted_percentages)
################################################################
file_path33 <- "./studies/deepCIS/ipmArthDAPmotifs/all_binding_sites.csv"
# Import the CSV file with all the bed-files TRUE TARGETS
true_binding_sites <- read.csv(file_path33)
################################################################
# Load the predicted percentages file
predicted_percentages <- read.csv("./studies/deepCIS/ipmArthDAPmotifs/predicted_percentages_per_family.csv")
###############################################################
# View the first few rows of the updated dataframe
head(true_binding_sites)
true_binding_sites$center <- c((((true_binding_sites$end+true_binding_sites$start)/2)+0.5))
true_binding_sites$wstart <- c(true_binding_sites$center-125)
true_binding_sites$wend <- c(true_binding_sites$center+125)
true_binding_sites$wlength <- c(true_binding_sites$wend-true_binding_sites$wstart) 
# Create the new column 'seq_id' by concatenating 'chr', 'start', and 'end' with the desired format
true_binding_sites$seq_id <- paste0(true_binding_sites$chr, ":", true_binding_sites$wstart, "-", true_binding_sites$wend)
###############################################################
common_ids <- read.table(file_path11, sep = " ", header = FALSE, stringsAsFactors = FALSE)
colnames(common_ids) <- c("seq_id", "label")
common_ids <- common_ids %>% distinct(seq_id,label, .keep_all = TRUE)
# Define the chunk size (number of rows to read at a time)
chunk_size <- 100000  # Adjust this number based on your memory capacity
# Count the total number of rows in the file
total_rows <- fread(file_path, select = 1)[, .N]
# Initialize an empty list to store chunks
chunks <- list()
chunk_number <- 1
skip <- 0
while(skip < total_rows) {
  # Calculate the number of rows to read
  rows_left <- total_rows - skip
  nrows_to_read <- min(chunk_size, rows_left)
  # Read a chunk of the file
  chunk <- fread(file_path, nrows = nrows_to_read, skip = skip, header = (chunk_number == 1))
  # Store the chunk in the list
  chunks[[chunk_number]] <- chunk
  # Increment the chunk number and skip value
  chunk_number <- chunk_number + 1
  skip <- skip + nrows_to_read
}
# Combine all chunks into one data.table (optional)
extracted_sites <- rbindlist(chunks)
# Display the first few rows of the combined data
head(extracted_sites)
# ####################################### WRANGLING PREPARE EXTR FILES
extr_seq_ids <- extracted_sites[, c("seq_id")]
# Split the `seq_id` column
extr_seq_ids <- extr_seq_ids %>%
  separate(seq_id, into = c("chr", "loc"), sep = ":") %>%
  separate(loc, into = c("start", "end"), sep = "-")

# Convert columns to numeric
extr_seq_ids <- extr_seq_ids %>%
  mutate(start = as.numeric(start),
         end = as.numeric(end))

# Add 1 to the `end` column
extr_seq_ids$end <- extr_seq_ids$end + 1
extr_seq_ids_reversed <- extr_seq_ids %>%
  mutate(seq_id = paste(chr, paste(start, end, sep = "-"), sep = ":"))

extracted_sites0 <- extracted_sites %>%
  select(-seq_id) %>%  # Remove the original `seq_id` column
  bind_cols(extr_seq_ids_reversed %>% select(seq_id))  # Add the new `seq_id` column


head(extracted_sites0)
#########################################
# Extract the first columns (IDs) from both dataframes
comm_ids_only <- common_ids$seq_id
extract_ids_only <- extracted_sites0$seq_id
head(comm_ids_only)
head(extract_ids_only)

# Find common IDs between both dataframes
true_ids <- intersect(extract_ids_only, comm_ids_only)
false_ids0 <- setdiff(extract_ids_only, comm_ids_only)
false_ids1 <- setdiff( comm_ids_only, extract_ids_only)
false_ids01 <- union(false_ids1, false_ids0)

#########################################
extracted_sites0_filtered <- extracted_sites0 %>% filter(seq_id %in% true_ids)
common_ids0_filtered <- common_ids %>% filter(seq_id %in% true_ids)

################################################################## This step generates the genome/chromosome wide binding and interatcion statistics
# Create a dataframe to store results
res_gen_wide_binding_stats <- data.frame(
  column_name = character(),
  total_extracted_sites = numeric(),
  total_common_ids = numeric(),
  predicted_binding = numeric(),
  motif_binding_context_value = numeric(),
  motif_redundancy_value = numeric(),
  true_binding_count = numeric(),        # New column for count of true binding
  true_per_total_extracted = numeric(),  # New column for intersection with extracted sites
  TDR_IPM = numeric(),          # New column for intersection with common ids
  FPFN_IPM_plus_binding_prediction = numeric(),
  stringsAsFactors = FALSE
)

# Loop over each column except 'seq_id'
for (col_name in colnames(extracted_sites0_filtered)) {
  if (col_name != "seq_id") {
    # Filter matching_common_ids for the current factor
    matching_common_ids <- common_ids0_filtered %>%
      filter(seq_id %in% extracted_sites0_filtered$seq_id &
               grepl(paste0("*", gsub("_tnt", "", col_name), "*"), label))
    
    # Filter the relevant column from extracted_sites0_filtered
    extracted_sites_filtered <- extracted_sites0_filtered %>%
      select(seq_id, all_of(col_name)) %>%
      distinct()
    
    # Perform semi-join to filter matching seq_ids
    extracted_sites_doublefiltered <- extracted_sites_filtered %>%
      semi_join(matching_common_ids, by = "seq_id")
    
    # Calculate the total number of extracted sites (total_mybs_extracted_sites)
    total_extracted_sites <- nrow(extracted_sites_doublefiltered)
    
    # Calculate the total number of common ids (total_mybs_common_ids)
    total_common_ids <- nrow(matching_common_ids)
    
    # Count values greater than 0.5 in the current column (count_mybs_above_0_5)
    predicted_binding <- sum(extracted_sites_doublefiltered[[col_name]] > 0.5)
    
    # Count values greater than 0.5 in the current column (count_mybs_above_0_5)
    predicted_binding_common <- sum(matching_common_ids[[col_name]] > 0.5)
    
    # Calculate the motif binding context value
    motif_binding_context_value <- 1-(predicted_binding/ total_extracted_sites)  
    
    # Calculate the ratio of total_common_ids to total_extracted_sites
    motif_redundancy_value <- ifelse(total_extracted_sites > 0, total_common_ids / total_extracted_sites, NA)
    
    # Make true_binding_count unique by using unique() to filter out duplicate seq_id values
    true_binding_count <- sum(
      unique(extracted_sites0_filtered$seq_id[extracted_sites0_filtered[[col_name]] > 0.5]) %in% true_binding_sites$seq_id
    )
    
    # Calculate true_per_total_extracted: intersection of seq_id between true_binding_sites and extracted_sites_doublefiltered
    true_per_total_extracted <- sum(extracted_sites_doublefiltered$seq_id %in% true_binding_sites$seq_id)
    
    TDR_IPM <- true_per_total_extracted/total_extracted_sites
    
    TDR_IPM_plus_binding_prediction <- true_binding_count/predicted_binding
    
    # Add the results to the dataframe
    res_gen_wide_binding_stats <- rbind(res_gen_wide_binding_stats, data.frame(
      column_name = col_name,
      total_extracted_sites = total_extracted_sites,
      total_common_ids = total_common_ids,
      predicted_binding = predicted_binding,
      motif_binding_context_value = motif_binding_context_value,
      motif_redundancy_value = motif_redundancy_value,
      true_binding_count = true_binding_count,
      true_per_total_extracted = true_per_total_extracted,
      TDR_IPM = TDR_IPM,
      FPFN_IPM_plus_binding_prediction = TDR_IPM_plus_binding_prediction
    ))
  }
}

# Print the results
head(res_gen_wide_binding_stats)

###################################################################
# Step 2: Rename the 'TF' column in the CSV to match 'column_name' in res_gen_wide_binding_stats
head(predicted_percentages)
# Forcefully clean column names of any unwanted characters
# Rename the first column directly using names()
names(predicted_percentages)[1] <- "column_name"
# Step 3: Merge the two dataframes based on 'column_name'
merged_stats <- merge(res_gen_wide_binding_stats, predicted_percentages, by = "column_name", all.x = TRUE)
# Step 4: Print the merged dataframe
print(merged_stats)
write.csv(merged_stats, file = output_file_path, row.names = FALSE)
##################################################################
# Assuming merged_stats is the dataframe
# Remove rows with NA values in 'motif_binding_context_value' and 'PredictedPercentage'
merged_stats_cleaned <- merged_stats %>%
  filter(!is.na(motif_binding_context_value) & !is.na(Predicted.Percentage))

# Fit a polynomial model (quadratic)
fit_poly <- lm(motif_binding_context_value ~ poly(Predicted.Percentage, 2), data = merged_stats_cleaned)

# Get the R-squared value
r_squared <- summary(fit_poly)$r.squared
print(paste("R-squared: ", r_squared))

# Generate predictions using the model
merged_stats_cleaned$predicted_values <- predict(fit_poly, newdata = data.frame(Predicted.Percentage = merged_stats_cleaned$Predicted.Percentage))

# Calculate the residuals (difference between observed and predicted values)
merged_stats_cleaned$residuals <- abs(merged_stats_cleaned$motif_binding_context_value - merged_stats_cleaned$predicted_values)

# Define outliers: points with residuals greater than a chosen threshold (e.g., 0.5 standard deviations from the mean residual)
threshold <- mean(merged_stats_cleaned$residuals) + 0.5 * sd(merged_stats_cleaned$residuals)
merged_stats_cleaned$outlier <- merged_stats_cleaned$residuals > threshold

# Calculate the standard deviation of 'motif_binding_context_value'
std_dev <- sd(merged_stats_cleaned$motif_binding_context_value)

# Create the ggplot
p <- ggplot(merged_stats_cleaned, aes(x = Predicted.Percentage, y = (1 / motif_binding_context_value)-1)) +
  # Scatter plot with transparency adjustment based on PredictedPercentage
  geom_point(aes(color = outlier, alpha = ifelse(Predicted.Percentage < 0.6, 0.4, 1)), size = 3) +  
  # Polynomial fit line
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, color = "darkgrey") +  
  # Label all points with transparency adjustment
  geom_text_repel(aes(label = str_remove(column_name, "_tnt$"), alpha = ifelse(Predicted.Percentage < 0.6, 0.4, 1)),
                  size = 4,                 # Increase label size
                  box.padding = 0.6,        # Increase padding around labels
                  point.padding = 0.3,      # Increase padding around points
                  max.overlaps = 15) +      # Allow more overlaps before labels are removed
  scale_alpha_identity() +  # Use the alpha values as set directly in the data
  scale_color_manual(values = c("black", "orange"), guide = FALSE) +  # Color outliers in orange
  labs(title = "",
       x = "Sensitivity",
       y = "Predictability (1/IPMciv)-1") +
  annotate("text", x = min(merged_stats_cleaned$Predicted.Percentage), y = max(1 / merged_stats_cleaned$motif_binding_context_value)-1,
           label = paste("R-squared:", round(r_squared, 3)), color = "darkgrey", hjust = 0) +  # R-squared in plot
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14, face = "bold"), # Larger axis labels
    axis.text = element_text(size = 12),                # Larger axis text
    axis.line = element_line(size = 0.8),               # Add thicker axis lines
    axis.ticks = element_line(size = 0.8),              # Add thicker axis ticks
    plot.title = element_text(size = 16, face = "bold") # Larger plot title
  )

# Print the plot in a new window
print(p)

# Save the plot as PNG
ggsave(filename = "./motif_binding_context_vs_predicted_percentage.png", plot = p, width =14, height = 6, dpi = 900)



##################################################################################################################
head(extracted_sites0_filtered)
data <- extracted_sites0_filtered
# Filter seq_ids for values greater than 0.5 per column
seq_id_lists <- lapply(data %>% select(-seq_id), function(col) {
  data$seq_id[which(col > 0.5)]
})

# Find overlapping seq_ids between columns
###############################
###############################
###############################
###############################
overlap_count <- matrix(0, ncol = length(seq_id_lists), nrow = length(seq_id_lists),
                        dimnames = list(names(seq_id_lists), names(seq_id_lists)))

for (i in 1:length(seq_id_lists)) {
  for (j in 1:length(seq_id_lists)) {
    if (i != j) {
      overlap_count[i, j] <- length(intersect(seq_id_lists[[i]], seq_id_lists[[j]]))
    }
  }
}
###############################
###############################
###############################
###############################
# Display the overlap count matrix
print(overlap_count)
# overlap_count
shared_seq_id <- as.data.frame(overlap_count)
colnames(overlap_count)
#### START HERE
shared_seq_id <- shared_seq_id %>%
  select(HB_tnt, SBP_tnt, BZR_tnt, LOBAS2_tnt, BBRBPC_tnt, AP2EREBP_tnt)  # Use column names without quotes
###############################
###############################
# Z-Score Normalization function
z_score_normalize <- function(x) {
  return((x - mean(x)) / sd(x))
}
###############################
###############################
# Apply Z-Score normalization to the original shared_seq_id
shared_seq_id_z_normalized <- shared_seq_id %>%
  mutate(across(everything(), z_score_normalize))
print(shared_seq_id_z_normalized)

# Filter rows where at least one value is above 1
filtered_shared_seq_id <- shared_seq_id_z_normalized %>%
  filter(apply(., 1, function(x) any(x > 1)))

# Modify column names to remove "_tnt"
colnames(filtered_shared_seq_id) <- str_remove(colnames(filtered_shared_seq_id), "_tnt$")
# Modify column names to remove "_tnt"
row.names(filtered_shared_seq_id) <- str_remove(row.names(filtered_shared_seq_id), "_tnt$")
# Check the filtered data
print(filtered_shared_seq_id)

# Convert to matrix for clustering
shared_seq_matrix <- as.matrix(filtered_shared_seq_id)

# Calculate distance and perform hierarchical clustering
dist_matrix <- dist(shared_seq_matrix, method = "euclidean")
hc_rows <- hclust(dist_matrix, method = "complete")  # Hierarchical clustering on rows
hc_cols <- hclust(dist(t(shared_seq_matrix)), method = "complete")  # Hierarchical clustering on columns

# Generate the heatmap
h <- pheatmap(shared_seq_matrix,
              cluster_rows = hc_rows,      # Use the clustering result for rows
              cluster_cols = hc_cols,      # Use the clustering result for columns
              display_numbers = TRUE,          # Display values on the heatmap
              fontsize_number = 8,             # Adjust number size
              main = "",
              fontsize = 10,                   # Adjust main font size
              color = colorRampPalette(c("white", "darkorange"))(100))  # Color gradient

print(h)

# Save the heatmap to a file
ggsave(filename = "./low_motif_binding_outliers_heatmap.png", plot = h, width = 6, height = 6, dpi = 900)
################################################################################
##################################################################################
#Build global network
# Display the overlap count matrix
print(overlap_count)
# overlap_count
shared_ALL_seq_id <- as.data.frame(overlap_count)
colnames(overlap_count)
###############################

# Define a safe Z-score normalization function
z_score_normalize_safe <- function(x) {
  if (sd(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x))) # Replace with zeros if standard deviation is zero
  } else {
    return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  }
}
###############################
###############################
###############################

# Clean the column names (remove any extra spaces)
predicted_percentages$TF <- gsub(" ", "", predicted_percentages$TF)

# Assuming shared_seq_id is already created in previous parts of the code.
# Add the predicted percentage values to the shared sequence IDs
shared_seq_id_with_percentage <- shared_seq_id %>%
  mutate(TF = rownames(shared_seq_id)) %>%
  left_join(predicted_percentages, by = "TF")

# Apply Z-Score normalization
shared_ALL_seq_id_id_z_normalized <- shared_ALL_seq_id %>%
  mutate(across(everything(), z_score_normalize_safe))

# Replace values greater than 1 with the original value, others with NA
filtered_data <- shared_ALL_seq_id_id_z_normalized %>%
  mutate(across(everything(), ~ ifelse(. > 1, ., NA)))

# Flatten the data for edge creation
edges <- filtered_data %>%
  as.data.frame() %>%
  rownames_to_column(var = "from") %>%
  pivot_longer(cols = -from, names_to = "to", values_to = "weight") %>%
  filter(!is.na(weight))

# Create a list of all unique nodes (TFs) from the edges and all available TFs
all_nodes <- unique(c(edges$from, edges$to))

# Create the graph object, ensuring all nodes are included, even those without edges
graph <- graph_from_data_frame(edges, vertices = data.frame(name = all_nodes))

# Check if all nodes are present in shared_seq_id_with_percentage
missing_nodes <- setdiff(V(graph)$name, shared_seq_id_with_percentage$TF)
if (length(missing_nodes) > 0) {
  message("Missing nodes in predicted percentages: ", paste(missing_nodes, collapse = ", "))
}

# Merge the predicted percentage data into the graph's vertices
# Ensure that the `Predicted_Percentage` values are assigned to matching nodes in the graph
V(graph)$Predicted_Percentage <- shared_seq_id_with_percentage$Predicted.Percentage[match(V(graph)$name, shared_seq_id_with_percentage$TF)]

# Check the length of V(graph)$Predicted_Percentage
if (length(V(graph)$Predicted_Percentage) != length(V(graph))) {
  stop("Mismatch in length of Predicted_Percentage and number of vertices in the graph.")
}

# Map Predicted_Percentage to node color using a gradient
V(graph)$node_color <- scales::col_numeric(palette = c("red", "blue"), domain = c(0, 1))(V(graph)$Predicted_Percentage)

# Remove the "_tnt" suffix from the labels
V(graph)$label <- str_remove(V(graph)$name, "_tnt$")
###########
# Compute degree (number of connections) for each node
degree_values <- degree(graph)

# Set the node size based on degree
node_size <- scales::rescale(degree_values, to = c(2, 10))  # Map degree to node size range

###########
# Plot the graph with the new coloring logic for nodes and light grey edges
N <- ggraph(graph, layout = "fr") +
  geom_edge_link(aes(width = weight), color = "darkgrey", alpha = 0.3) +  # Set edges to light grey
  geom_node_point(aes(color = Predicted_Percentage, size = node_size)) +  # Map Predicted_Percentage for color
  #and node_size for node size
  geom_node_text(aes(label = label), repel = TRUE, size = 4, box.padding = 0.4) +  # Adjust text for clarity
  theme_void() +
  scale_color_gradient(low = "red", high = "blue", name = "Sensitivity") +  # Continuous color gradient from red to blue
  scale_size_continuous(name = "Centrality", range = c(2, 10)) +  # Map node size to degree, range from 3 to 10
  ggtitle("") +
  guides(color = guide_colorbar(title = "Sensitivity", barwidth = 1, barheight = 6),  # Adjust color legend width
         size = guide_legend(title = "Centrality", override.aes = list(color = "black"))) # Add legend for node size

# Print the plot
print(N)

# Save the plot to a file
ggsave(filename = "Filtered_Zscore_Network_Optimized_with_legend.png", plot = N, width = 7, height = 9, dpi = 900)

# Create a gradient from red to blue with purple in the middle
color_palette <- colorRampPalette(c("red", "purple", "blue"))

# Generate the gradient image
gradient <- matrix(seq(0, 85, length.out = 256), nrow = 1)

# Plot the gradient
image(gradient, col = color_palette(256), axes = FALSE, xlab = "", ylab = "")
# Add color bar
colorbar <- colorRampPalette(c("red", "purple", "blue"))(256)
axis(1, at = c(0, 1), labels = c("0", "85"))


###############################
###############################
###############################
###############################
# Apply Z-Score normalization
shared_ALL_seq_id_id_z_normalized <- shared_ALL_seq_id %>%
  mutate(across(everything(), z_score_normalize_safe))

# Replace values greater than 1 with the original value, others with NA
filtered_data <- shared_ALL_seq_id_id_z_normalized %>%
  mutate(across(everything(), ~ ifelse(. > 1, ., NA)))


head(filtered_data)
# Define `node_names` based on unique names in the filtered data
# Define `node_names` based on the column names in `filtered_data`
node_names <- colnames(filtered_data)

# Trim `_tnt` suffix from `node_names`
node_names <- str_remove(node_names, "_tnt$")

# Define central nodes and trim `_tnt`
central_node <- c("BZR_tnt", "BES1_tnt", "bHLH_tnt", "bZIP_tnt", "SBP_tnt") %>%
  str_remove("_tnt$")

# Create edges connecting only to central nodes with Z > 1 and assign edge colors
edges <- bind_rows(
  data.frame(
    from = "BZR_tnt",
    to = node_names[!is.na(filtered_data$BZR_tnt)],
    weight = filtered_data$BZR_tnt[!is.na(filtered_data$BZR_tnt)]
  ),
  data.frame(
    from = "BES1_tnt",
    to = node_names[!is.na(filtered_data$BES1_tnt)],
    weight = filtered_data$BES1_tnt[!is.na(filtered_data$BES1_tnt)]
  ),
  data.frame(
    from = "bHLH_tnt",
    to = node_names[!is.na(filtered_data$bHLH_tnt)],
    weight = filtered_data$bHLH_tnt[!is.na(filtered_data$bHLH_tnt)]
  ),
  data.frame(
    from = "bZIP_tnt",
    to = node_names[!is.na(filtered_data$bZIP_tnt)],
    weight = filtered_data$bZIP_tnt[!is.na(filtered_data$bZIP_tnt)]
  ),
  data.frame(
    from = "SBP_tnt",
    to = node_names[!is.na(filtered_data$SBP_tnt)],
    weight = filtered_data$SBP_tnt[!is.na(filtered_data$SBP_tnt)]
  )
)

# Modify edge and node labels by removing the "_tnt" suffix
edges <- edges %>%
  mutate(
    from = str_remove(from, "_tnt$"),
    to = str_remove(to, "_tnt$")
  )

# Assign colors to edges based on connection type
edges <- edges %>%
  mutate(
    color = case_when(
      from == "BZR" | to == "BZR" ~ "orange",        # Edges connected to BZR
      from == "bHLH" | to == "bHLH" ~ "aquamarine", # Edges connected to bHLH
      from == "BES1" | to == "BES1" ~ "plum",       # Edges connected to BES1
      TRUE ~ "grey"                                 # Other edges
    )
  )

# Ensure node names are unique and valid
all_node_names <- unique(c(central_node, edges$to))

# Create the graph object
graph <- graph_from_data_frame(edges, vertices = data.frame(name = all_node_names), directed = FALSE)

# Plot the graph
N <- ggraph(graph, layout = "fr") +
  geom_edge_link(aes(width = weight, color = color), alpha = 0.6) +
  geom_node_point(size = 5, aes(color = case_when(
    name == "bHLH" ~ "aquamarine3", 
    name == "BES1" ~ "plum3", 
    name %in% c("SBP", "BZR") ~ "darkorange", 
    TRUE ~ "black"
  ))) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_edge_color_manual(values = c("orange" = "orange", "aquamarine" = "aquamarine", "plum" = "plum", "grey" = "grey")) + # Apply edge colors
  scale_color_identity() + # Use colors directly assigned
  theme_void() +
  theme(legend.title = element_blank()) +
  ggtitle("")
print(N)

# Save the plot
ggsave(filename = "G-Box_TF_Network_Relevant_Nodes_Full_Network.png", plot = N, width = 4, height = 3, dpi = 900)



###############################

###############################
###############################
#rm(list = ls())  # Clear all objects from memory
#gc()             # Force garbage collection
# Apply Z-score normalization
shared_ALL_seq_id_id_z_normalized <- shared_ALL_seq_id %>%
  mutate(across(everything(), z_score_normalize_safe))

# Filter data for Z-scores > 1
filtered_data <- shared_ALL_seq_id_id_z_normalized %>%
  mutate(across(everything(), ~ ifelse(. > 1, ., NA)))

# Ensure no values <= 1 remain
print(filtered_data)

# Flatten into a long-format table for edges
edges <- filtered_data %>%
  as.data.frame() %>%
  rownames_to_column(var = "from") %>%
  pivot_longer(cols = -from, names_to = "to", values_to = "weight") %>%
  filter(!is.na(weight))  # Retain only edges with Z > 1

# Check if invalid edges exist
invalid_edges <- edges %>% filter(from == "HB_tnt" & to == "CPP")
print(invalid_edges)  # Should show 0 rows
# Debug edges
print(edges %>% filter(weight <= 1))  # Ensure no invalid edges
# Clear old graph objects if necessary
rm(graph)
# Modify edge and node labels by removing the "_tnt" suffix
edges <- edges %>%
  mutate(
    from = str_remove(from, "_tnt$"),
    to = str_remove(to, "_tnt$")
  )

# Update node names in the graph object to remove "_tnt"
all_node_names <- unique(c(edges$from, edges$to))
graph <- graph_from_data_frame(
  d = edges,
  vertices = data.frame(name = str_remove(all_node_names, "_tnt$")),
  directed = FALSE
)
# Assign colors to edges
relevant_nodes <- c("BZR", "SBP", "HB")
# Reapply edge colors safely
edges <- edges %>%
  mutate(
    color = ifelse(from %in% relevant_nodes | to %in% relevant_nodes, "orange", "grey")
  )
tail(edges)
# Ensure edge width reflects correct weight values
geom_edge_link(aes(width = weight, color = color), alpha = 0.6)
# Apply a minimum width threshold to edges

# Check for unexpected edge colors or NA
print(edges %>% filter(is.na(color)))

# Create graph object
all_node_names <- unique(c(edges$from, edges$to))
# Recreate the graph object
graph <- graph_from_data_frame(edges, vertices = data.frame(name = unique(c(edges$from, edges$to))), directed = FALSE)
#graph <- graph_from_data_frame(edges, vertices = data.frame(name = all_node_names), directed = FALSE)


M <- ggraph(graph, layout = "fr") +
  geom_edge_link(aes(width = weight, color = color), alpha = 0.6) +
  geom_node_point(size = 2, aes(color = ifelse(name %in% relevant_nodes, "highlighted", "other"))) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_edge_color_manual(
    name = "Edge Types",
    values = c("orange" = "orange", "grey" = "grey"),
    labels = c("orange" = "low IPMciv outliers", "grey" = "others")
  ) +
  scale_edge_width(
    name = "Z-scores",  # Rename weight legend to "Z-scores"
    range = c(0.5, 2)   # Adjust width scaling (optional)
  ) +
  scale_color_manual(
    name = "Node Types",
    values = c("highlighted" = "darkorange", "other" = "black"),
    labels = c("highlighted" = "low IPMciv outliers", "other" = "others")
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.spacing = unit(1, "lines"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  ) +
  ggtitle("")

print(M)
print(edges %>% filter(from == "HB_tnt" & to == "CPP"))

# Save the plot
ggsave(filename = "Filtered_Zscore_Network_Relevant_Nodes_Full_Network.png", plot = M, width = 6, height = 8, dpi = 900)
#################################################################################
#################################################################################

###############################
# Apply Z-Score normalization
shared_ALL_seq_id_id_z_normalized <- shared_ALL_seq_id %>%
  mutate(across(everything(), z_score_normalize_safe))

# Replace values greater than 1 with the original value, others with NA
filtered_data <- shared_ALL_seq_id_id_z_normalized %>%
  mutate(across(everything(), ~ ifelse(. > 1, ., NA)))

# Define `node_names` based on the column names in `filtered_data`
node_names <- colnames(filtered_data)

# Trim `_tnt` suffix from `node_names`
node_names <- str_remove(node_names, "_tnt$")

# Define central nodes and trim `_tnt`
central_node <- c("ARID_tnt", "Homeobox_tnt", "Homeobox_ecoli", "ZFHD_tnt", "MYBrelated_tnt", "HB_tnt", "CPP_tnt", "C2C2dof_tnt", "C2C2YABBY_tnt", "REM_tnt") %>%
  str_remove("_tnt$")

# Create edges connecting only to central nodes with Z > 1 and assign edge colors
edges <- bind_rows(
  data.frame(
    from = "ARID_tnt",
    to = node_names[!is.na(filtered_data$ARID_tnt)],
    weight = filtered_data$ARID_tnt[!is.na(filtered_data$ARID_tnt)]
  ),
  data.frame(
    from = "Homeobox_tnt",
    to = node_names[!is.na(filtered_data$Homeobox_tnt)],
    weight = filtered_data$Homeobox_tnt[!is.na(filtered_data$Homeobox_tnt)]
  ),
  data.frame(
    from = "Homeobox_ecoli",
    to = node_names[!is.na(filtered_data$Homeobox_ecoli)],
    weight = filtered_data$Homeobox_ecoli[!is.na(filtered_data$Homeobox_ecoli)]
  ),
  data.frame(
    from = "ZFHD_tnt",
    to = node_names[!is.na(filtered_data$ZFHD_tnt)],
    weight = filtered_data$ZFHD_tnt[!is.na(filtered_data$ZFHD_tnt)]
  ),
  data.frame(
    from = "MYBrelated_tnt",
    to = node_names[!is.na(filtered_data$MYBrelated_tnt)],
    weight = filtered_data$MYBrelated_tnt[!is.na(filtered_data$MYBrelated_tnt)]
  ),
  data.frame(
    from = "HB_tnt",
    to = node_names[!is.na(filtered_data$HB_tnt)],
    weight = filtered_data$HB_tnt[!is.na(filtered_data$HB_tnt)]
  ),
  data.frame(
    from = "CPP_tnt",
    to = node_names[!is.na(filtered_data$CPP_tnt)],
    weight = filtered_data$CPP_tnt[!is.na(filtered_data$CPP_tnt)]
  ),
  data.frame(
    from = "C2C2dof_tnt",
    to = node_names[!is.na(filtered_data$C2C2dof_tnt)],
    weight = filtered_data$C2C2dof_tnt[!is.na(filtered_data$C2C2dof_tnt)]
  ),
  data.frame(
    from = "C2C2YABBY_tnt",
    to = node_names[!is.na(filtered_data$C2C2YABBY_tnt)],
    weight = filtered_data$C2C2YABBY_tnt[!is.na(filtered_data$C2C2YABBY_tnt)]
  ),
  data.frame(
    from = "REM_tnt",
    to = node_names[!is.na(filtered_data$REM_tnt)],
    weight = filtered_data$REM_tnt[!is.na(filtered_data$REM_tnt)]
  )
)

# Modify edge and node labels by removing the "_tnt" suffix
edges <- edges %>%
  mutate(
    from = str_remove(from, "_tnt$"),
    to = str_remove(to, "_tnt$")
  )

# Assign colors to edges based on connection type
edges <- edges %>%
  mutate(
    color = ifelse(from == "HB" | to == "HB", "orange", "grey")
  )

# Ensure node names are unique and valid
all_node_names <- unique(c(central_node, edges$to))

# Create the graph object
graph <- graph_from_data_frame(edges, vertices = data.frame(name = all_node_names), directed = FALSE)

# Plot the graph
N <- ggraph(graph, layout = "fr") +
  geom_edge_link(aes(width = weight, color = color), alpha = 0.6) +
  geom_node_point(size = 5, color = "black") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_edge_color_manual(
    name = "Edge Types",
    values = c("orange" = "orange", "grey" = "grey"),
    labels = c("orange" = "Connections to HB", "grey" = "Other Connections")
  ) +
  theme_void() +
  ggtitle("")
print(N)

# Save the plot
ggsave(filename = "Subnetwork_Specified_Factors_HB_Highlighted.png", plot = N, width = 6, height = 6, dpi = 900)
