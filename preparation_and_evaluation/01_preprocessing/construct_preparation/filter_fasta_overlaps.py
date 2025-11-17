import argparse
import re
import sys

def parse_header(header_line):
    """
    Parses a FASTA header to extract chromosome and start position.
    Expected format: >Chr:Start-End ...
    Example: >1:10029917-10030166 ...
    Returns: (chromosome, start_position) or (None, None) if parsing fails.
    """
    match = re.match(r">\s*([^:]+):(\d+)-(\d+)", header_line)
    if match:
        chromosome = match.group(1)
        start_pos = int(match.group(2))
        # We primarily care about the start position for the overlap check
        return chromosome, start_pos
    else:
        print(f"Warning: Could not parse header format: {header_line.strip()}", file=sys.stderr)
        return None, None

def filter_overlapping_sequences(input_fasta, output_fasta, threshold):
    """
    Reads sequences from input_fasta, filters based on coordinate overlap,
    and writes the kept sequences to output_fasta.

    Args:
        input_fasta (str): Path to the input FASTA file.
        output_fasta (str): Path to the output FASTA file.
        threshold (int): The +/- range to consider sequences overlapping based
                         on their start coordinates.
    """
    kept_regions = {}  # Dictionary to store {chromosome: [list of kept start positions]}
    current_header = None
    current_sequence = []
    sequences_processed = 0
    sequences_kept = 0

    print(f"Processing {input_fasta}...")
    print(f"Overlap threshold: +/- {threshold} bp from start position.")

    try:
        with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
            for line in infile:
                line = line.strip()
                if not line:  # Skip empty lines
                    continue

                if line.startswith(">"):
                    sequences_processed += 1
                    # --- Process the *previous* sequence before starting the new one ---
                    if current_header:
                        chrom, start = parse_header(current_header)
                        if chrom is not None:
                            is_overlapping = False
                            # Check for overlap only if chromosome has been seen before
                            if chrom in kept_regions:
                                for kept_start in kept_regions[chrom]:
                                    if abs(start - kept_start) <= threshold:
                                        is_overlapping = True
                                        print(f"  Skipping: {current_header.split()[0]} (overlaps with region starting near {kept_start} on chr {chrom})", file=sys.stderr)
                                        break # Found an overlap, no need to check further

                            if not is_overlapping:
                                # Keep this sequence
                                sequences_kept += 1
                                outfile.write(current_header + "\n")
                                outfile.write("\n".join(current_sequence) + "\n")
                                # Add its start position to the list of kept regions
                                if chrom not in kept_regions:
                                    kept_regions[chrom] = []
                                kept_regions[chrom].append(start)
                        else:
                             # Could not parse header, decide if you want to keep it anyway
                             # Defaulting to skipping sequences with unparseable headers
                             print(f"  Skipping sequence with unparseable header: {current_header}", file=sys.stderr)


                    # --- Start the new sequence ---
                    current_header = line
                    current_sequence = []
                elif current_header: # Only append if we are currently inside a valid sequence entry
                    current_sequence.append(line)

            # --- Process the *last* sequence in the file ---
            if current_header:
                chrom, start = parse_header(current_header)
                if chrom is not None:
                    is_overlapping = False
                    if chrom in kept_regions:
                        for kept_start in kept_regions[chrom]:
                            if abs(start - kept_start) <= threshold:
                                is_overlapping = True
                                print(f"  Skipping: {current_header.split()[0]} (overlaps with region starting near {kept_start} on chr {chrom})", file=sys.stderr)
                                break

                    if not is_overlapping:
                        sequences_kept += 1
                        outfile.write(current_header + "\n")
                        outfile.write("\n".join(current_sequence) + "\n")
                        # Optionally add to kept_regions if needed for further checks,
                        # but not strictly necessary after the last sequence
                        # if chrom not in kept_regions:
                        #     kept_regions[chrom] = []
                        # kept_regions[chrom].append(start)
                else:
                    print(f"  Skipping last sequence with unparseable header: {current_header}", file=sys.stderr)


    except FileNotFoundError:
        print(f"Error: Input file not found: {input_fasta}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"\nProcessing complete.")
    print(f"Total sequences read: {sequences_processed}")
    print(f"Sequences kept (non-overlapping): {sequences_kept}")
    print(f"Output written to: {output_fasta}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter a FASTA file to remove sequences whose start coordinates "
                    "overlap within a given threshold (+/- bp) of a previously kept sequence. "
                    "Keeps the first sequence encountered in an overlapping set."
    )
    parser.add_argument("input_fasta", help="Path to the input FASTA file.")
    parser.add_argument("output_fasta", help="Path to the output filtered FASTA file.")
    parser.add_argument(
        "-t", "--threshold", type=int, default=40,
        help="The +/- threshold distance between start coordinates to define overlap (default: 40)."
    )

    args = parser.parse_args()

    if args.threshold < 0:
        print("Error: Threshold must be a non-negative integer.", file=sys.stderr)
        sys.exit(1)

    filter_overlapping_sequences(args.input_fasta, args.output_fasta, args.threshold)
