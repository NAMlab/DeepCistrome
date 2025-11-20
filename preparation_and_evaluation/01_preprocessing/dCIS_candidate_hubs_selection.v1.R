################################################################################
# install.packages("wordcloud")
library(wordcloud)
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)

################################################################################
setwd("/home/ibg-4/Desktop/deepCIS_calc/")
file_path <- "./studies/deepCIS/ipmArthDAPmotifs/extracted_sites_predictions_chrom_1.csv"
file_path11 <- "./studies/deepCIS/ipmArthDAPmotifs/extracted_ranges_1.txt"
################################################################################
selection <- c("TCP",
               "WRKY",
               "bHLH",
               "BES1",
               "bZIP",
               "Trikelix",
               "MYB",
               "CAMPTRA",
               "FAR1",
               "MYB")
################################################################################
################################################################################

common_ids <- read.table(file_path11, sep = " ", header = FALSE, stringsAsFactors = FALSE)
colnames(common_ids) <- c("seq_id", "label")
common_ids <- common_ids %>% distinct(seq_id,label, .keep_all = TRUE)

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
extracted_ids0_filtered_2 <- extracted_sites0_filtered[rowSums(extracted_sites0_filtered > 0.5) >= 1, ] # FILTER FOR SITES THAT ARE PREDICTED TO BE BOUND BY TWO FACTORS MIN
# colnames(common_ids0_filtered_2)
# Identify columns that match any string in the selection list
selected_columns <- grep(paste(selection, collapse = "|"), names(extracted_ids0_filtered_2), value = TRUE)

# Filter rows where no value in the selected columns is above 0.5
extracted_sites0_filtered_3 <- extracted_ids0_filtered_2[rowSums(extracted_ids0_filtered_2[, ..selected_columns] > 0.5) > 0, ]
################################################################################ # FILTER FOR SITES THAT ARE NOT INTERSECTING

################################################################################
head(common_ids0_filtered)
head(extracted_sites0_filtered_3)
merge00 <- merge(common_ids0_filtered, extracted_sites0_filtered_3, by = "seq_id")
merge00 <- unique(merge00)
merge00$count_above_05 <- rowSums(merge00[, selected_columns] > 0.5)
print(summary(merge00$count_above_05))
print(unique(merge00$count_above_05))

#################################################################################
#################################################################################
merge01 <- merge00 %>%
  rowwise() %>%
  filter(any(sapply(selection, function(s) grepl(s, label, fixed = TRUE)))) %>%
  ungroup()


# Extract the substring between the second and third underscores
merge01 <- merge01 %>%
  mutate(target_column = str_match(label, "^[^_]+_[^_]+_([^_]+)_")[, 2])

colnames(merge01)
head(merge01$target_column)
# Filter rows where the value in the target column is greater than 0.5
################################################################################ Keep only windows where IPM and PREDCTION MACTCHES. 
merge02 <- merge01 %>%
  rowwise() %>%
  filter(any(sapply(names(.), function(colname) {
    grepl(target_column, colname, fixed = TRUE) && get(colname) > 0.5
  }))) %>%
  ungroup()
################################################################################
# Extract chromosome, start, and end positions from seq_id
merge03 <- merge02 %>%
  separate(seq_id, into = c("chr", "start", "end"), sep = ":|-", convert = TRUE)

# Sort the data by chromosome and start position
merge03 <- merge03 %>%
  arrange(chr, start)

# Function to check for overlaps using a sweep line algorithm
find_non_overlapping <- function(data) {
  non_overlapping <- c()
  current_end <- -Inf
  
  for (i in seq_len(nrow(data))) {
    if (data$start[i] > current_end || data$chr[i] != data$chr[i - 1]) {
      non_overlapping <- c(non_overlapping, i)
      current_end <- data$end[i]
    }
  }
  
  return(non_overlapping)
}
# Get indices of non-overlapping ranges
non_overlapping_indices <- find_non_overlapping(merge03)

# Filter the data to keep only non-overlapping ranges
merge04 <- merge03[non_overlapping_indices, ]
merge04 <- merge04 %>%
  mutate(seq_id = paste(chr, start, sep = ":")) %>%
  mutate(seq_id = paste0(seq_id, "-", end-1))

head(merge04$seq_id)
################################################################################
# Filter rows where count_above_05 equals 6
bin_6_data <- merge04[merge04$count_above_05 == 6 & 5 & 4 & 3 & 2, ]

# Select the relevant columns for the tag cloud
selected_columns <- grep(paste(selection, collapse = "|"), names(bin_6_data), value = TRUE)

# Create a word frequency table
word_freq <- sort(colSums(bin_6_data[, selected_columns, drop = FALSE]), decreasing = TRUE)

# Generate the tag cloud
wordcloud(names(word_freq), word_freq, scale = c(4, 0.5),
          min.freq = 1, max.words = 200, random.order = FALSE, rot.per = 0.35,
          colors = brewer.pal(8, "Dark2"))
# Plot the histogram with correct breaks
# Calculate the frequency of each count (excluding 0)
count_table <- table(merge04$count_above_05)

# Remove the count for 0 if it exists
count_table <- count_table[count_table > 0]

# Create a bar plot
bar_positions <- barplot(count_table,
                         main = "Distribution of Binding Hubs",
                         xlab = "Number of Predicted Binding Partners",
                         ylab = "Frequency",
                         col = "skyblue",
                         border = "black",
                         names.arg = names(count_table))  # Ensure correct bin labeling

text(x = barplot(count_table, plot = FALSE),
     y = count_table,
     labels = count_table,
     pos = 1,  # Position labels at the top of the bars
     cex = 1.2,
     col = "red",
     offset = -1)  # Adjust the offset to ensure labels are visible
################################################################################
head(merge02$seq_id)
# Save the DataFrame as a CSV file
write.csv(merge04, "dCIS_candidate_hubs_background.csv", row.names = FALSE)

# Export the seq_id column as a list and save it as a text file
seq_id_list <- merge04$seq_id
write.table(seq_id_list, "dCIS_candidate_hubs_seqID_list.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

