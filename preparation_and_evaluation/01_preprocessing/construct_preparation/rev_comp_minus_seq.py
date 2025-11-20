#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import textwrap # For wrapping sequence lines

def reverse_complement(seq):
    """
    Computes the reverse complement of a DNA sequence.

    Handles standard bases (ATCG), ambiguity codes (NRYWSMK),
    and preserves case.

    Args:
        seq (str): The DNA sequence.

    Returns:
        str: The reverse complement of the sequence.
    """
    # Define complement mapping (includes ambiguity codes & case)
    # R(AG) Y(CT) W(AT) S(GC) M(AC) K(GT)
    complement_map = str.maketrans(
        "ATCGNRYWSMKatcgnrywsmk",
        "TAGCNYRWSMKtagcnyrwsmk"
    )
    # Reverse the sequence string
    reversed_seq = seq[::-1]
    # Apply the complement map
    return reversed_seq.translate(complement_map)

def process_fasta_for_revcomp(input_filepath, output_stream):
    """
    Reads a fasta file, reverse complements sequences with '(-)' in the header,
    and writes to the output stream.

    Args:
        input_filepath (str): Path to the input fasta file.
        output_stream (file object): Where to write the output (e.g., sys.stdout or a file).
    """
    header = None
    sequence_parts = []
    line_count = 0
    entries_processed = 0
    entries_revcomped = 0

    try:
        with open(input_filepath, 'r') as infile:
            for line in infile:
                line_count += 1
                line = line.strip()
                if not line:
                    continue # Skip empty lines

                if line.startswith('>'):
                    # Process the previous entry if one exists
                    if header:
                        entries_processed += 1
                        sequence = "".join(sequence_parts)
                        # Check if header contains the minus orientation marker
                        if "(-)" in header:
                            new_sequence = reverse_complement(sequence)
                            entries_revcomped += 1
                        else:
                            new_sequence = sequence

                        # Write the header
                        print(header, file=output_stream)
                        # Write the (potentially modified) sequence, wrapped
                        # Use textwrap for cleaner wrapping
                        wrapped_sequence = textwrap.fill(new_sequence, width=60)
                        print(wrapped_sequence, file=output_stream)

                        sequence_parts = [] # Reset for next entry

                    # Store the new header
                    header = line
                else:
                    # Append sequence line if we are currently processing an entry
                    if header:
                        sequence_parts.append(line)
                    else:
                        # Sequence line before any header - issue with file format?
                        print(f"Warning: Found sequence line before first header near line {line_count}. Skipping.", file=sys.stderr)


            # Process the very last entry in the file after the loop ends
            if header:
                entries_processed += 1
                sequence = "".join(sequence_parts)
                # Check if header contains the minus orientation marker
                if "(-)" in header:
                    new_sequence = reverse_complement(sequence)
                    entries_revcomped += 1
                else:
                    new_sequence = sequence

                # Write the header
                print(header, file=output_stream)
                # Write the (potentially modified) sequence, wrapped
                wrapped_sequence = textwrap.fill(new_sequence, width=60)
                print(wrapped_sequence, file=output_stream)

        print(f"\nProcessed {entries_processed} fasta entries.", file=sys.stderr)
        print(f"Reverse complemented {entries_revcomped} entries marked with '(-)' in header.", file=sys.stderr)

    except FileNotFoundError:
        print(f"Error: Input file not found at {input_filepath}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)


# --- Main Execution ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Reverse complements sequences in a multifasta file if their "
                    "header line contains the marker '(-)' (e.g., from 'iur(-)' or 'igr(-)' tags).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input_fasta",
                        help="Path to the input multifasta file.")
    parser.add_argument("-o", "--output",
                        help="Optional path to write the output fasta file. "
                             "If not specified, output goes to standard output (console).")

    args = parser.parse_args()

    # Determine output stream
    output_stream = sys.stdout
    output_file_handle = None
    if args.output:
        try:
            output_file_handle = open(args.output, 'w')
            output_stream = output_file_handle
            print(f"Output will be written to: {args.output}", file=sys.stderr)
        except IOError as e:
            print(f"Error: Could not open output file {args.output} for writing: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        print("Output will be written to standard output.", file=sys.stderr)

    # Process the fasta file
    process_fasta_for_revcomp(args.input_fasta, output_stream)

    # Close the output file if it was opened
    if output_file_handle:
        output_file_handle.close()
        print(f"Output successfully written to {args.output}", file=sys.stderr)

    print("\nScript finished.", file=sys.stderr)
