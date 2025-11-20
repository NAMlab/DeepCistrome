import argparse
import re
import sys
from collections import defaultdict

# --- Configuration ---
# Define the substrings that uniquely identify control and variant headers.
# Adjust these if your actual headers use different identifiers.
CTRL_IDENTIFIER = "dCIS-ctrl"
VAR_IDENTIFIER = "dCIS-var"
# --- End Configuration ---

def parse_coord_string(header_line):
    """
    Parses a FASTA header to extract the coordinate string "Chr:Start-End".
    Returns: The coordinate string or None if parsing fails.
    """
    match = re.match(r">\s*([^:]+:\d+-\d+)", header_line)
    if match:
        # Validate start <= end within the coordinate string itself
        parts = re.match(r"([^:]+):(\d+)-(\d+)", match.group(1))
        if parts:
            try:
                start = int(parts.group(2))
                end = int(parts.group(3))
                if start <= end:
                    return match.group(1) # Return the full "chr:start-end"
                else:
                    print(f"Warning: Header has start > end, skipping: {header_line.strip()}", file=sys.stderr)
                    return None
            except ValueError:
                 print(f"Warning: Non-integer coordinates in header: {header_line.strip()}", file=sys.stderr)
                 return None
        else: # Should not happen if outer regex matched, but belt-and-suspenders
             print(f"Warning: Could not re-parse coordinate string: {match.group(1)}", file=sys.stderr)
             return None

    else:
        print(f"Warning: Could not parse coordinate string from header: {header_line.strip()}", file=sys.stderr)
        return None

def get_entry_type(header_line):
    """
    Determines if the header represents a 'ctrl' or 'var' entry based on identifiers.
    """
    if CTRL_IDENTIFIER in header_line:
        return 'ctrl'
    elif VAR_IDENTIFIER in header_line:
        return 'var'
    else:
        return 'unknown'

def filter_paired_sequences(input_fasta, output_fasta):
    """
    Reads input_fasta, keeps only sequences that form a ctrl/var pair
    based on identical coordinate strings in their headers, and writes
    the paired sequences to output_fasta.
    """
    # Dictionary to store info grouped by coordinate string
    # Format: {'chr:start-end': {'has_ctrl': bool, 'has_var': bool, 'entries': [...]}}
    # where 'entries' stores {'header': str, 'sequence': str} dicts
    coords_info = defaultdict(lambda: {'has_ctrl': False, 'has_var': False, 'entries': []})

    current_header = None
    current_sequence_parts = []
    entries_read = 0

    # --- Pass 1: Read file, parse info, and group by coordinates ---
    print("Pass 1: Reading input file and grouping by coordinates...")
    try:
        with open(input_fasta, 'r') as infile:
            for line in infile:
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    # Process previous entry
                    if current_header is not None:
                        entries_read += 1
                        coord_string = parse_coord_string(current_header)
                        entry_type = get_entry_type(current_header)
                        sequence = "".join(current_sequence_parts)

                        if coord_string and sequence and entry_type != 'unknown':
                            # Store entry info, grouped by coordinate string
                            info = coords_info[coord_string]
                            info['entries'].append({'header': current_header, 'sequence': sequence})
                            if entry_type == 'ctrl':
                                info['has_ctrl'] = True
                            elif entry_type == 'var':
                                info['has_var'] = True
                        elif not coord_string:
                            print(f"Info: Skipping entry due to unparsed coordinate: {current_header}", file=sys.stderr)
                        elif not sequence:
                             print(f"Info: Skipping entry with empty sequence: {current_header}", file=sys.stderr)
                        elif entry_type == 'unknown':
                             print(f"Info: Skipping entry with unknown type (no '{CTRL_IDENTIFIER}' or '{VAR_IDENTIFIER}'): {current_header}", file=sys.stderr)


                    # Start new entry
                    current_header = line
                    current_sequence_parts = []
                elif current_header:
                    current_sequence_parts.append(line)

            # Process the last entry
            if current_header is not None:
                entries_read += 1
                coord_string = parse_coord_string(current_header)
                entry_type = get_entry_type(current_header)
                sequence = "".join(current_sequence_parts)

                if coord_string and sequence and entry_type != 'unknown':
                    info = coords_info[coord_string]
                    info['entries'].append({'header': current_header, 'sequence': sequence})
                    if entry_type == 'ctrl':
                        info['has_ctrl'] = True
                    elif entry_type == 'var':
                        info['has_var'] = True
                # Print warnings/info messages similar to above if needed for the last entry

    except FileNotFoundError:
        print(f"Error: Input file not found: '{input_fasta}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred during Pass 1: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Pass 1 complete. Read {entries_read} entries, grouped into {len(coords_info)} unique coordinate ranges.")

    # --- Pass 2: Filter coordinate groups and write paired sequences ---
    print(f"Pass 2: Writing entries that have both '{CTRL_IDENTIFIER}' and '{VAR_IDENTIFIER}' types for the same coordinates...")
    sequences_written = 0
    pairs_found = 0
    try:
        with open(output_fasta, 'w') as outfile:
            for coord_string, info in coords_info.items():
                # Check if this coordinate group has BOTH ctrl and var types
                if info['has_ctrl'] and info['has_var']:
                    pairs_found += 1
                    # Write all entries associated with this coordinate string
                    for entry in info['entries']:
                        outfile.write(entry['header'] + "\n")
                        outfile.write(entry['sequence'] + "\n")
                        sequences_written += 1
                # else: # Optional: Log discarded coordinates
                    # print(f"  Discarding coordinate group (missing ctrl/var pair): {coord_string}", file=sys.stderr)

    except Exception as e:
        print(f"An error occurred during Pass 2 (writing output): {e}", file=sys.stderr)
        sys.exit(1)

    print("\nProcessing complete.")
    print(f"Coordinate ranges with pairs found: {pairs_found}")
    print(f"Total sequences written (from paired groups): {sequences_written}")
    print(f"Output written to '{output_fasta}'")


# --- Main execution block ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=f"Filter a multi-FASTA file, keeping only sequences that form a\n"
                    f"control/variant pair based on identical coordinate strings \n"
                    f"(e.g., >Chr:Start-End) in their headers and the presence of \n"
                    f"'{CTRL_IDENTIFIER}' and '{VAR_IDENTIFIER}' identifiers."
    )
    parser.add_argument("input_fasta", help="Path to the input multi-FASTA file.")
    parser.add_argument("output_fasta", help="Path for the output filtered FASTA file containing only paired sequences.")

    args = parser.parse_args()
    filter_paired_sequences(args.input_fasta, args.output_fasta)
