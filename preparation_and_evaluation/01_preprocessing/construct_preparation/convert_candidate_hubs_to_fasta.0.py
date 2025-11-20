import pandas as pd

# Define file paths
insertion_file_path = "insertion_site_sequence.txt"
# Assuming the CSV file is named TFBS_mutagenesis-2.csv based on your 'head' command
csv_file_path = "TFBS_mutagenesis-2.csv"
output_fasta_path = "TFBS_mutated_dCIS_candidates.fasta"

# --- Step 1: Read insertion site sequence and split into flanks ---
try:
    with open(insertion_file_path, 'r') as f:
        insertion_template = f.read().strip()
except FileNotFoundError:
    print(f"Error: Insertion site file not found at '{insertion_file_path}'")
    exit()

placeholder = "<core promoter insertion>"
if placeholder not in insertion_template:
    print(f"Error: Placeholder '{placeholder}' not found in '{insertion_file_path}'")
    exit()

parts = insertion_template.split(placeholder)
flank1 = parts[0]
flank2 = parts[1]

print(f"Flank1 loaded (length: {len(flank1)}): {flank1[:30]}...")
print(f"Flank2 loaded (length: {len(flank2)}): {flank2[:30]}...")

# --- Step 2: Read CSV and process each row ---
try:
    df = pd.read_csv(csv_file_path)
    print(f"\nSuccessfully read '{csv_file_path}'. Processing {len(df)} rows.")
except FileNotFoundError:
    print(f"Error: CSV file not found at '{csv_file_path}'")
    exit()
except Exception as e:
    print(f"Error reading CSV file '{csv_file_path}': {e}")
    exit()

# Open output FASTA file
with open(output_fasta_path, 'w') as fasta_out:
    processed_count = 0
    # Iterate through each row in the DataFrame
    for index, row in df.iterrows():
        # --- 2.1: Build FASTA header ---
        try:
            sp = row['sp']
            gene = row['gene']
            chromosome = row['chr'] # Column name is 'chr' in the CSV head
            start_pos = row['start']
            end_pos = row['end']
            strand = row['strand']
            motif = row['motif']
            variant = row['variant']
            csv_sequence = str(row['sequence']) # Ensure sequence is a string
        except KeyError as e:
            print(f"\nError: Missing expected column in CSV: {e} at row {index+1}.")
            print(f"Please check your CSV file headers. Available columns: {df.columns.tolist()}")
            print("Skipping this row.")
            continue
        except TypeError as e:
            print(f"\nError: Type error for data in row {index+1}, likely in 'sequence' column (expected string): {e}")
            print("Skipping this row.")
            continue


        # Format the header as per specifications
        fasta_header = f">{sp}_{gene}_chr{chromosome}_{start_pos}..{end_pos}_({strand})_{motif}_{variant}"

        # --- 2.2: Build FASTA sequence ---
        # Insert CSV sequence into flanks
        full_constructed_seq = flank1 + csv_sequence + flank2

        # Trim the sequence to 250 bp from the center
        target_len = 250
        current_len = len(full_constructed_seq)

        final_sequence = ""
        if current_len < target_len:
            print(f"Warning: Row {index+1} (Gene: {gene}, Motif: {motif}): Constructed sequence length ({current_len}bp) is less than target {target_len}bp. Sequence will not be trimmed and will be shorter than {target_len}bp.")
            final_sequence = full_constructed_seq
        else:
            trim_total = current_len - target_len
            trim_from_start = trim_total // 2
            # Slice the sequence to get the central 250bp
            final_sequence = full_constructed_seq[trim_from_start : trim_from_start + target_len]

        # Write to FASTA file
        fasta_out.write(fasta_header + "\n")
        fasta_out.write(final_sequence + "\n")
        processed_count += 1

print(f"\nProcessed {processed_count} entries.")
print(f"Multi-FASTA file '{output_fasta_path}' created successfully.")
