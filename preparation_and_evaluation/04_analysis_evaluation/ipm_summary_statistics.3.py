# Import necessary libraries
import pandas as pd
import numpy as np

# Set the working directory (optional)
import os
os.chdir("/home/ibg-4/Desktop/Rhome/deepCIS_calc/workflows/")

# Print starting message
print("Starting script...")

# Define the file paths
file_path = "./../studies/deepCIS/ipmArthDAPmotifs/extracted_sites_predictions_chrom_1.csv"
file_path11 = "./../studies/deepCIS/ipmArthDAPmotifs/extracted_ranges_1.txt"
file_path22 = "./../studies/deepCIS/ipmArthDAPmotifs/predicetd_percentages_per_family.csv"
file_path33 = "./../studies/deepCIS/ipmArthDAPmotifs/all_binding_sites.csv"

# Load the CSV with predicted percentages
print("Loading predicted percentages...")
predicted_percentages = pd.read_csv(file_path22)

# Load the true binding sites
print("Loading true binding sites...")
true_binding_sites = pd.read_csv(file_path33)

# Calculate the 'center', 'wstart', 'wend', and 'wlength' columns
print("Calculating center, wstart, and wend for true binding sites...")
true_binding_sites['center'] = ((true_binding_sites['start'] + true_binding_sites['end']) / 2) + 0.5
true_binding_sites['wstart'] = true_binding_sites['center'] - 125
true_binding_sites['wend'] = true_binding_sites['center'] + 125
true_binding_sites['seq_id'] = true_binding_sites['chr'].astype(str) + ":" + true_binding_sites['wstart'].astype(int).astype(str) + "-" + true_binding_sites['wend'].astype(int).astype(str)

# Load the common_ids file
print("Loading common IDs...")
common_ids = pd.read_csv(file_path11, sep=" ", header=None, names=["seq_id", "label"])
common_ids = common_ids.drop_duplicates(subset=["seq_id", "label"])

# Load the extracted_sites_predictions_chrom_1 file in chunks
print("Loading extracted sites predictions...")
extracted_sites = pd.read_csv(file_path)

# Split the 'seq_id' into 'chr', 'start', 'end'
print("Splitting 'seq_id' into 'chr', 'start', and 'end'...")
extracted_sites[['chr', 'start_end']] = extracted_sites['seq_id'].str.split(':', expand=True)
extracted_sites[['start', 'end']] = extracted_sites['start_end'].str.split('-', expand=True)

# Convert 'start' and 'end' to numeric
print("Converting 'start' and 'end' to numeric and adjusting...")
extracted_sites['start'] = pd.to_numeric(extracted_sites['start'])
extracted_sites['end'] = pd.to_numeric(extracted_sites['end']) + 1

# Recreate 'seq_id' based on adjusted 'start' and 'end'
extracted_sites['seq_id'] = extracted_sites['chr'] + ":" + extracted_sites['start'].astype(str) + "-" + extracted_sites['end'].astype(str)

# Filter true_ids and false_ids
print("Filtering true and false IDs...")
comm_ids_only = common_ids['seq_id']
extract_ids_only = extracted_sites['seq_id']

true_ids = np.intersect1d(extract_ids_only, comm_ids_only)

# Filter extracted_sites0_filtered and common_ids0_filtered based on true_ids
print("Filtering extracted sites and common IDs based on true IDs...")
extracted_sites0_filtered = extracted_sites[extracted_sites['seq_id'].isin(true_ids)]
common_ids0_filtered = common_ids[common_ids['seq_id'].isin(true_ids)]

# Initialize the dataframe for storing results
print("Initializing result storage...")
res_gen_wide_binding_stats = []

# Extract numeric columns (the columns of interest for calculations)
numeric_columns = extracted_sites0_filtered.select_dtypes(include=[np.number]).columns

# Loop through each numeric column
print("Starting calculations for each column...")
for i, col_name in enumerate(numeric_columns):
    print(f"Processing column {i + 1}/{len(numeric_columns)}: {col_name}")
    
    # Filter matching_common_ids for the current factor
    factor_name = col_name.replace('_tnt', '')
    matching_common_ids = common_ids0_filtered[common_ids0_filtered['label'].str.contains(factor_name, regex=True)]
    
    # Filter the relevant column from extracted_sites0_filtered
    extracted_sites_filtered = extracted_sites0_filtered[['seq_id', col_name]].drop_duplicates()
    
    # Perform semi-join to filter matching seq_ids
    extracted_sites_doublefiltered = pd.merge(extracted_sites_filtered, matching_common_ids, on='seq_id', how='inner')

    # Calculate statistics
    total_extracted_sites = len(extracted_sites_doublefiltered)
    total_common_ids = len(matching_common_ids)
    count_above_0_5 = (extracted_sites_doublefiltered[col_name] > 0.5).sum()
    motif_binding_context_value = count_above_0_5 / total_extracted_sites if total_extracted_sites > 0 else np.nan
    motif_redundancy_value = total_common_ids / total_extracted_sites if total_extracted_sites > 0 else np.nan
    
    # Calculate true_binding_count
    true_binding_count = sum(extracted_sites_doublefiltered.loc[extracted_sites_doublefiltered[col_name] > 0.5, 'seq_id'].isin(true_binding_sites['seq_id']))
    
    # Calculate true_per_total_extracted
    true_per_total_extracted = sum(extracted_sites_doublefiltered['seq_id'].isin(true_binding_sites['seq_id']))
    
    # Calculate true_per_commons
    true_per_commons = sum(matching_common_ids['seq_id'].isin(true_binding_sites['seq_id']))
    
    # Append the result to the list
    res_gen_wide_binding_stats.append({
        'column_name': col_name,
        'total_extracted_sites': total_extracted_sites,
        'total_common_ids': total_common_ids,
        'count_above_0_5': count_above_0_5,
        'motif_binding_context_value': motif_binding_context_value,
        'motif_redundancy_value': motif_redundancy_value,
        'true_binding_count': true_binding_count,
        'true_per_total_extracted': true_per_total_extracted,
        'true_per_commons': true_per_commons
    })

print("Finished calculations. Converting results to DataFrame...")

# Convert the list to a dataframe
res_gen_wide_binding_stats = pd.DataFrame(res_gen_wide_binding_stats)

# Merge with predicted_percentages
print("Merging results with predicted percentages...")
predicted_percentages.rename(columns={'TF': 'column_name'}, inplace=True)
merged_stats = pd.merge(res_gen_wide_binding_stats, predicted_percentages, on='column_name', how='left')

# Save the merged_stats dataframe to a CSV file
output_file_path = "/home/ibg-4/Desktop/Rhome/deepCIS_calc/workflows/merged_stats_results.csv"
print(f"Saving results to {output_file_path}...")
merged_stats.to_csv(output_file_path, index=False)

# Print the final merged dataframe (optional)
print("Script completed. Results saved.")
