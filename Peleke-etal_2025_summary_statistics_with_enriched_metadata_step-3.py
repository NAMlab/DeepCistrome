import pandas as pd
import os
import re
import argparse
import sys
# Process TFBS predictions and enrich with JASPAR metadata.
#    Enrich the motif map track (.bed-file) summary.tsv with metadata using a database .transfac file. 
#    "Input-1: Large TFBS predictions .tsv file"
#    "Input-2: JASPAR .transfac file"
#    "Output 1: Summary statistics file"
#    "Output 2: Enriched statistics file"
# ==========================================
# PART 1: PROCESS LARGE FILE & GENERATE SUMMARY
# ==========================================
def process_predictions_in_chunks(input_file, summary_output_file, chunk_size=100000):
    """
    Reads large TSV in chunks, aggregates counts, saves the summary file,
    and returns the DataFrame for further processing.
    """
    print(f"[INFO] Starting chunked processing of {input_file}...")
    
    total_ones = None
    total_rows = None
    
    # Check delimiter
    sep = '\t' if input_file.endswith('.tsv') else ','
    
    try:
        reader = pd.read_csv(input_file, sep=sep, chunksize=chunk_size)
    except Exception as e:
        sys.exit(f"[ERROR] Could not read file: {e}")

    for i, chunk in enumerate(reader):
        print(f" > Processing chunk {i+1}...", end='\r')
        
        # 1. Parse TF_id
        if 'SequenceHeader' in chunk.columns:
            chunk['TF_id'] = chunk['SequenceHeader'].astype(str).apply(lambda x: x.split('_')[0])
        else:
            chunk['TF_id'] = chunk.iloc[:, 0].astype(str).apply(lambda x: x.split('_')[0])
            
        # 2. Identify numeric columns (DeepCIS Families)
        numeric_cols = chunk.select_dtypes(include=['number']).columns
        
        # 3. Aggregate
        chunk_sum = chunk.groupby('TF_id')[numeric_cols].sum()   # Sum of 1s
        chunk_count = chunk.groupby('TF_id')[numeric_cols].count() # Total rows
        
        # 4. Update global totals
        if total_ones is None:
            total_ones = chunk_sum
            total_rows = chunk_count
        else:
            total_ones = total_ones.add(chunk_sum, fill_value=0)
            total_rows = total_rows.add(chunk_count, fill_value=0)

    print(f"\n[INFO] Aggregation complete. Calculating statistics...")

    # Fill NaNs
    total_ones = total_ones.fillna(0)
    total_rows = total_rows.fillna(0)
    
    # Calculate Zeros and Recovery Rate
    total_zeros = total_rows - total_ones
    recovery_rates = total_ones / total_rows

    # 5. Reshape (Stacking)
    s_ones = total_ones.stack().rename('Count_1')
    s_zeros = total_zeros.stack().rename('Count_0')
    s_total = total_rows.stack().rename('Total')
    # Renamed from 'Ratio' to 'Recovery_Rate'
    s_rate = recovery_rates.stack().rename('Recovery_Rate')

    # Combine
    summary_df = pd.concat([s_zeros, s_ones, s_total, s_rate], axis=1)
    summary_df.index.names = ['TF_id', 'Predicted_Family']

    # 6. Save Output File #1 (The Summary)
    print(f"[INFO] Saving summary statistics to {summary_output_file}...")
    summary_df.to_csv(summary_output_file, sep='\t')
    
    return summary_df

# ==========================================
# PART 2: METADATA PARSING
# ==========================================
def parse_jaspar_metadata(file_path):
    """
    Parses JASPAR TRANSFAC file for metadata.
    """
    metadata = {}
    current_id = None
    current_data = {}
    
    patterns = {
        'tf_name': re.compile(r'^ID\s+(.*)'), 
        'JASPAR_tf_family': re.compile(r'CC\s+tf_family:(.*)'),
        'tf_class': re.compile(r'CC\s+tf_class:(.*)'),
        'data_type': re.compile(r'CC\s+data_type:(.*)'),
        'pubmed_ids': re.compile(r'CC\s+pubmed_ids:(.*)'),
        'uniprot_ids': re.compile(r'CC\s+uniprot_ids:(.*)')
    }

    print(f"[INFO] Parsing metadata from {file_path}...")
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('AC'):
                if current_id and current_data:
                    metadata[current_id] = current_data
                parts = line.split()
                current_id = parts[1].strip() if len(parts) > 1 else None
                current_data = {k: 'N/A' for k in patterns.keys()}
            
            if current_id:
                for key, pattern in patterns.items():
                    if key == 'tf_name' and line.startswith('ID'):
                        match = pattern.match(line)
                        if match: current_data[key] = match.group(1).strip()
                    elif line.startswith('CC'):
                        match = pattern.match(line)
                        if match: current_data[key] = match.group(1).split(';')[0].strip()
            
            if line.startswith('//'):
                if current_id and current_data:
                    metadata[current_id] = current_data
                    current_id = None
                    current_data = {}
                    
        if current_id and current_data:
            metadata[current_id] = current_data

    return metadata

# ==========================================
# PART 3: MAIN PIPELINE
# ==========================================
def main():
    parser = argparse.ArgumentParser(description="Process TFBS predictions and enrich with JASPAR metadata.")
    parser.add_argument("-p", "--predictions", required=True, help="Input: Large TFBS predictions .tsv file")
    parser.add_argument("-m", "--metadata", required=True, help="Input: JASPAR .transfac file")
    parser.add_argument("--summary_out", required=True, help="Output 1: Summary statistics file")
    parser.add_argument("--enriched_out", required=True, help="Output 2: Enriched statistics file")
    parser.add_argument("--split_dir", default="split_results", help="Output Directory: Folder to store the split class/family files")
    
    args = parser.parse_args()

    # 1. Validation
    if not os.path.exists(args.predictions):
        sys.exit(f"Error: Prediction file '{args.predictions}' not found.")
    if not os.path.exists(args.metadata):
        sys.exit(f"Error: Metadata file '{args.metadata}' not found.")

    # 2. Process Large File -> Save Summary -> Get DataFrame
    summary_df = process_predictions_in_chunks(args.predictions, args.summary_out)

    # 3. Prepare DataFrame for Enrichment
    # Reset index to make 'TF_id' and 'Predicted_Family' normal columns
    df_for_enrichment = summary_df.reset_index()

    # 4. Load & Parse Metadata
    meta_dict = parse_jaspar_metadata(args.metadata)
    meta_df = pd.DataFrame.from_dict(meta_dict, orient='index')
    meta_df.index.name = 'TF_id'
    
    # Prepare Base ID for fallback (MA0001.1 -> MA0001)
    df_for_enrichment['Base_ID'] = df_for_enrichment['TF_id'].astype(str).apply(lambda x: x.split('.')[0])
    meta_df['Base_ID'] = meta_df.index.astype(str).map(lambda x: x.split('.')[0])

    # 5. Merge
    print("[INFO] Merging statistics with metadata...")
    merged = df_for_enrichment.merge(meta_df.drop(columns=['Base_ID']), on='TF_id', how='left')

    # 6. Fill Missing Metadata (Base ID Fallback)
    missing_mask = merged['JASPAR_tf_family'].isna() | (merged['JASPAR_tf_family'] == 'N/A')
    if missing_mask.sum() > 0:
        print(f"[INFO] {missing_mask.sum()} rows missing metadata. Using Base ID fallback...")
        meta_base_lookup = meta_df.reset_index().drop_duplicates(subset='Base_ID').set_index('Base_ID')
        if 'TF_id' in meta_base_lookup.columns:
            meta_base_lookup = meta_base_lookup.drop(columns=['TF_id'])
            
        fallback_df = df_for_enrichment.merge(meta_base_lookup, on='Base_ID', how='left')
        
        cols = ['tf_name', 'JASPAR_tf_family', 'tf_class', 'pubmed_ids', 'uniprot_ids', 'data_type']
        for col in cols:
            if col in merged.columns and col in fallback_df.columns:
                merged[col] = merged[col].fillna(fallback_df[col])

    # 7. Cleanup and Order
    if 'Base_ID' in merged.columns:
        merged = merged.drop(columns=['Base_ID'])
        
    # --- REORDERING COLUMNS ---
    # Section 1: Metadata (Left)
    meta_cols = ['TF_id', 'tf_name', 'JASPAR_tf_family', 'tf_class', 'data_type', 'pubmed_ids', 'uniprot_ids']
    # Section 2: DeepCIS Predictions (Middle)
    prediction_cols = ['Predicted_Family']
    # Section 3: Statistics (Right)
    stat_cols = ['Count_1', 'Count_0', 'Total', 'Recovery_Rate']
    
    # Combine
    desired_order = meta_cols + prediction_cols + stat_cols
    
    # Select only columns that actually exist (to handle potential missing metadata fields)
    final_order = [c for c in desired_order if c in merged.columns]
    
    # Add any remaining columns that weren't in our explicit list at the end (safeguard)
    remaining = [c for c in merged.columns if c not in final_order]
    final_order = final_order + remaining
    
    merged = merged[final_order]

    # 8. Save Output File #2 (Enriched)
    print(f"[INFO] Saving enriched statistics to {args.enriched_out}...")
    merged.to_csv(args.enriched_out, sep='\t', index=False)

    # ==========================================
    # PART 9: GENERATE SEPARATE SUMMARY FILES
    # ==========================================
    print(f"[INFO] Splitting results by TF_class and Predicted_Family into '{args.split_dir}'...")
    
    # Ensure output directory exists
    if not os.path.exists(args.split_dir):
        os.makedirs(args.split_dir)

    # Handle NaNs before grouping so they don't get dropped
    merged['tf_class'] = merged['tf_class'].fillna('Unclassified')
    merged['Predicted_Family'] = merged['Predicted_Family'].fillna('Unknown')

    # Group by the required columns
    grouped = merged.groupby(['tf_class', 'Predicted_Family'])

    count_saved = 0
    for (tf_cls, pred_fam), group_df in grouped:
        # Sanitize filename (remove / or spaces)
        # e.g., "AP2/EREBP" becomes "AP2_EREBP"
        safe_cls = str(tf_cls).replace('/', '_').replace(' ', '_')
        safe_fam = str(pred_fam).replace('/', '_').replace(' ', '_')
        
        filename = f"{safe_cls}_{safe_fam}.tsv"
        filepath = os.path.join(args.split_dir, filename)
        
        group_df.to_csv(filepath, sep='\t', index=False)
        count_saved += 1
        
    print(f"[INFO] Done! Saved {count_saved} separate summary files in '{args.split_dir}'.")

if __name__ == "__main__":
    main()
