#!/usr/bin/env python3

import argparse
import sys
import os

def read_fasta(filepath):
    """
    Reads a multi-FASTA file and yields records as (header, sequence) tuples.
    Handles multi-line sequences.
    """
    header = None
    sequence_parts = []
    try:
        with open(filepath, 'r') as infile:
            for line in infile:
                line = line.strip()
                if not line: # Skip empty lines
                    continue
                if line.startswith('>'):
                    # If we have a previous record, yield it
                    if header is not None:
                        yield header, "".join(sequence_parts)
                    # Start new record
                    header = line[1:] # Store header without '>'
                    sequence_parts = []
                elif header is not None: # Only append if we are inside a record
                    # Append sequence line (remove potential whitespace within sequence?)
                    # For now, assume sequence lines have no internal whitespace
                    sequence_parts.append(line)
            # Yield the last record in the file
            if header is not None:
                yield header, "".join(sequence_parts)
    except FileNotFoundError:
        print(f"Error: Input FASTA file not found at '{filepath}'", file=sys.stderr)
        raise # Re-raise the exception to be caught in main
    except Exception as e:
        print(f"Error reading FASTA file '{filepath}': {e}", file=sys.stderr)
        raise

def write_fasta_record(outfile, header, sequence, wrap=70):
    """Writes a single FASTA record to the output file, wrapping sequence."""
    outfile.write(f">{header}\n")
    for i in range(0, len(sequence), wrap):
        outfile.write(sequence[i:i+wrap] + "\n")

def process_fasta(input_fasta_path, output_fasta_path):
    """
    Reads FASTA, filters, trims, sorts, and writes processed sequences.

    Args:
        input_fasta_path (str): Path to the input multi-FASTA file.
        output_fasta_path (str): Path for the output processed multi-FASTA file.

    Returns:
        bool: True if successful, False otherwise.
    """
    # --- Strings to filter sequences by ---
    filter_strings = ["GGTCTC", "GAGACC"]
    # --- Trimming parameters ---
    target_length = 170
    original_expected_length = 250
    trim_start_index = 40 # 0-based index to start keeping (removes first 40)
    trim_end_index = trim_start_index + target_length # Index to stop keeping (keeps up to index 209) -> 210

    print(f"Input FASTA: {input_fasta_path}")
    print(f"Output FASTA: {output_fasta_path}")
    print(f"Filtering out sequences containing: {', '.join(filter_strings)}")
    print(f"Trimming sequences to {target_length}bp (keeping base {trim_start_index + 1} through {trim_end_index})")

    initial_records = []
    filtered_records = []
    trimmed_records = []
    read_count = 0
    filter_count = 0
    trim_skip_count = 0
    trim_count = 0

    # --- Step 1: Read all records ---
    print("\nReading input FASTA...")
    try:
        initial_records = list(read_fasta(input_fasta_path))
        read_count = len(initial_records)
        print(f"Read {read_count} sequences.")
        if read_count == 0:
             print("Warning: No sequences found in input file.")
             # Write an empty output file? Or exit? Let's exit cleanly.
             # Create empty file to signify completion with no output
             open(output_fasta_path, 'w').close()
             return True

    except Exception as e:
        # Error already printed in read_fasta
        return False

    # --- Step 2: Filter sequences containing forbidden strings ---
    print("\nFiltering sequences...")
    for header, sequence in initial_records:
        found_filter_string = False
        for forbidden in filter_strings:
            if forbidden in sequence:
                print(f"  Filtering out sequence '{header}' (contains '{forbidden}')")
                filter_count += 1
                found_filter_string = True
                break # No need to check other forbidden strings for this sequence
        if not found_filter_string:
            filtered_records.append((header, sequence))

    print(f"Kept {len(filtered_records)} sequences after filtering (removed {filter_count}).")

    if not filtered_records:
        print("Warning: No sequences remained after filtering.")
        open(output_fasta_path, 'w').close()
        return True

    # --- Step 3: Trim sequences ---
    print("\nTrimming sequences...")
    for header, sequence in filtered_records:
        original_len = len(sequence)
        # Check if sequence is long enough to be trimmed correctly
        if original_len >= trim_end_index: # Need at least 210 bases to trim correctly
             if original_len != original_expected_length:
                  print(f"  Warning: Sequence '{header}' has length {original_len} (expected {original_expected_length}), but trimming anyway.")
             trimmed_sequence = sequence[trim_start_index:trim_end_index]
             trimmed_records.append((header, trimmed_sequence))
             trim_count += 1
             # Verify trimmed length (optional check)
             # if len(trimmed_sequence) != target_length:
             #    print(f"  Warning: Trimming '{header}' resulted in length {len(trimmed_sequence)}, expected {target_length}.")
        else:
            print(f"  Warning: Sequence '{header}' has length {original_len}, which is too short to trim to {target_length}bp (needs at least {trim_end_index}bp). Skipping.")
            trim_skip_count += 1

    print(f"Trimmed {trim_count} sequences to {target_length}bp (skipped {trim_skip_count} due to insufficient length).")

    if not trimmed_records:
        print("Warning: No sequences remained after trimming.")
        open(output_fasta_path, 'w').close()
        return True

    # --- Step 4: Sort by header ---
    print("\nSorting remaining sequences by header...")
    trimmed_records.sort(key=lambda item: item[0])
    print("Sorting complete.")

    # --- Step 5: Write output ---
    print(f"\nWriting {len(trimmed_records)} processed sequences to '{output_fasta_path}'...")
    try:
        with open(output_fasta_path, 'w') as outfile:
            for header, sequence in trimmed_records:
                write_fasta_record(outfile, header, sequence)
        print("Writing complete.")
        return True
    except IOError as e:
        print(f"Error writing to output file '{output_fasta_path}': {e}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"An unexpected error occurred during writing: {e}", file=sys.stderr)
        return False


# --- Main execution block ---
if __name__ == "__main__":
    # --- Argument Parser Setup ---
    parser = argparse.ArgumentParser(
         description="Processes a multi-FASTA file: filters out sequences containing specific strings, "
                     "trims sequences to a central 170bp region (from original 250bp), "
                     "sorts by header, and writes the result.",
         formatter_class=argparse.ArgumentDefaultsHelpFormatter
         )

    # --- Define Arguments ---
    parser.add_argument(
        "input_fasta",
        type=str,
        help="Path to the input multi-FASTA file (expecting sequences of approx. 250bp)."
    )
    parser.add_argument(
        "output_fasta",
        type=str,
        help="Path for the output processed multi-FASTA file."
    )

    # --- Parse the arguments ---
    args = parser.parse_args()

    # --- Basic Check ---
    # Input file existence check is handled within read_fasta

    # --- Run the main logic ---
    if not process_fasta(args.input_fasta, args.output_fasta):
         print("\nScript finished with errors.", file=sys.stderr)
         sys.exit(1) # Exit with error code
    else:
         print("\nScript finished successfully.")
         sys.exit(0) # Exit with success code
