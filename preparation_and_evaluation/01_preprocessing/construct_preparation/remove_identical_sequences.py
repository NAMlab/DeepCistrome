import argparse
import sys

def process_fasta_sequences(input_file, output_file):
    """
    Reads a multi-FASTA file, removes entries with duplicate sequences
    (keeping the first occurrence), and writes the unique entries to an output file.
    Sequence comparison is case-insensitive.

    Args:
        input_file (str): Path to the input multi-FASTA file.
        output_file (str): Path to the output multi-FASTA file.
    """
    seen_sequences = set()
    current_header = None
    current_sequence_parts = []
    sequences_read = 0
    sequences_written = 0

    print(f"Processing '{input_file}'...")

    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                line = line.strip()
                if not line: # Skip empty lines
                    continue

                if line.startswith(">"):
                    # Process the previously gathered sequence (if any)
                    if current_header is not None:
                        sequences_read += 1
                        # Join sequence parts and convert to uppercase for case-insensitive comparison
                        full_sequence = "".join(current_sequence_parts).upper()

                        # Check if the sequence is non-empty and hasn't been seen before
                        if full_sequence and full_sequence not in seen_sequences:
                            # Add to seen set
                            seen_sequences.add(full_sequence)
                            # Write the kept sequence (original header, joined sequence parts)
                            outfile.write(current_header + "\n")
                            # Re-join original sequence parts without case change for writing
                            outfile.write("".join(current_sequence_parts) + "\n")
                            sequences_written += 1
                        # else: # Optional: uncomment to see which headers are skipped
                            # print(f"  Skipping duplicate sequence under header: {current_header}", file=sys.stderr)


                    # Start the new sequence entry
                    current_header = line
                    current_sequence_parts = []

                # Only append sequence lines if we have encountered a header
                elif current_header is not None:
                    # Append the raw sequence line (no case change yet)
                    current_sequence_parts.append(line)

            # --- Process the very last sequence in the file ---
            if current_header is not None:
                sequences_read += 1
                full_sequence = "".join(current_sequence_parts).upper()
                if full_sequence and full_sequence not in seen_sequences:
                    seen_sequences.add(full_sequence)
                    outfile.write(current_header + "\n")
                    outfile.write("".join(current_sequence_parts) + "\n")
                    sequences_written += 1
                # else: # Optional: uncomment to see which headers are skipped
                    # print(f"  Skipping duplicate sequence under header: {current_header}", file=sys.stderr)


    except FileNotFoundError:
        print(f"Error: Input file not found: '{input_file}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred during processing: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"\nProcessing complete.")
    print(f"Total sequences read: {sequences_read}")
    print(f"Unique sequences written: {sequences_written}")
    print(f"Output written to: '{output_file}'")

# --- Main execution block ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Remove duplicate sequences from a multi-FASTA file based on "
                    "100% sequence identity (case-insensitive), keeping the first "
                    "occurrence of each unique sequence."
    )
    parser.add_argument("input_fasta", help="Path to the input multi-FASTA file.")
    parser.add_argument("output_fasta", help="Path where the filtered output FASTA file will be saved.")

    # Parse arguments from the command line
    args = parser.parse_args()

    # Run the processing function
    process_fasta_sequences(args.input_fasta, args.output_fasta)
