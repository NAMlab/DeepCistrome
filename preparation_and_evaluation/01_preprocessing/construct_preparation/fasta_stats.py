#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import collections
import re

# --- Function Definitions ---

def generate_statistics(input_filepath):
    """
    Parses a tagged fasta file and generates statistics per TFBS family.
    'single'/'multi' counts based on number found after tag block in header.

    Args:
        input_filepath (str): Path to the input fasta file.
                              Assumed header format: >coords tag_block number TFBS_Family [opt]

    Returns:
        dict: A nested dictionary with statistics.
              Outer keys: TFBS Family
              Inner keys: 'total', 'iur', 'igr', '(?)', 'plus', 'minus', 'single', 'multi'
              Returns None if file cannot be read or no stats generated.
    """
    # Define the categories for counting, matching the desired output table
    categories = ['total', 'iur', 'igr', '(?)', 'plus', 'minus', 'single', 'multi']
    # Initialize statistics dictionary using defaultdict for convenience
    stats = collections.defaultdict(lambda: {category: 0 for category in categories})

    # Regex to capture coordinates and the rest of the line
    # Capture group 1: coords (e.g., 1:100-200)
    # Capture group 2: everything after the first space following coords
    header_regex_simple = re.compile(r"^>(\S+:\d+-\d+)\s+(.*)")

    print(f"Processing file: {input_filepath}...", file=sys.stderr)
    processed_headers = 0
    skipped_headers = 0

    try:
        with open(input_filepath, 'r') as infile:
            for line_num, line in enumerate(infile, 1):
                if not line.startswith('>'):
                    continue # Skip sequence lines

                processed_headers += 1
                header_content = line[1:].strip() # Remove '>' and surrounding whitespace

                # --- Parse Header Components ---
                coords_part = "N/A"
                tags_part = "N/A"
                number_str = None # To store the number string (e.g., "1", "2")
                tfbs_family = "Unknown" # Default TFBS family

                match = header_regex_simple.match(line)
                if match:
                    coords_part = match.group(1)
                    # Split the rest of the line after coordinates into max 2 parts:
                    # part 1 = tag block, part 2 = everything after tag block
                    tags_and_following = match.group(2).split(maxsplit=1)

                    if len(tags_and_following) >= 1:
                        tags_part = tags_and_following[0] # First part is always the tag block

                        if len(tags_and_following) == 2: # If there's more after the tag block
                             # Split the rest by whitespace to get number, TFBS, etc.
                             following_parts = tags_and_following[1].split()
                             if len(following_parts) >= 1:
                                 # Assume first item is the number
                                 number_str = following_parts[0]
                                 if len(following_parts) >= 2:
                                     # Assume second item is TFBS family
                                     tfbs_family = following_parts[1]
                                 # If only number found, tfbs_family remains "Unknown"
                             # If following_parts is empty, number/TFBS remain defaults
                        # If only tag block found, number/TFBS remain defaults
                    else:
                         # Only coordinates found after '>', edge case
                         tags_part = "N/A"
                else:
                     # Header didn't match basic pattern
                     print(f"Warning: Skipping malformed header on line {line_num}: {line.strip()}", file=sys.stderr)
                     skipped_headers += 1
                     continue # Skip this header
                # --- End Header Parsing ---


                # --- Incrementing Statistics ---
                # Get the dictionary for the current TFBS family (creates if needed)
                current_stats = stats[tfbs_family]
                current_stats['total'] += 1

                # Increment counts based on content of the 'tags_part'
                if "iur" in tags_part:
                    current_stats['iur'] += 1
                if "igr" in tags_part:
                    current_stats['igr'] += 1
                if tags_part == "(?)": # Check for exact unknown tag
                   current_stats['(?)'] += 1
                if "(+)" in tags_part: # Check for presence of plus orientation
                    current_stats['plus'] += 1
                if "(-)" in tags_part: # Check for presence of minus orientation
                    current_stats['minus'] += 1

                # Increment 'single'/'multi' based on the extracted number string
                if number_str is not None:
                    try:
                        num = int(number_str) # Convert string to integer
                        if num == 1:
                            current_stats['single'] += 1
                        elif 2 <= num <= 6: # Check if number is between 2 and 6 inclusive
                            current_stats['multi'] += 1
                        # Numbers outside 1-6 are ignored for these columns
                    except ValueError:
                        # The string found wasn't a valid integer
                        print(f"Warning: Expected number after tag block, but found '{number_str}' on line {line_num}: {line.strip()}", file=sys.stderr)
                # else: number_str was None, couldn't count single/multi
                # --- End Statistics Incrementing ---

        print(f"Processed {processed_headers} headers, skipped {skipped_headers}.", file=sys.stderr)
        if not stats:
            print("Warning: No statistics generated. Check file format and content.", file=sys.stderr)
            return None
        return stats

    except FileNotFoundError:
        print(f"Error: Input file not found at {input_filepath}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred during processing: {e}", file=sys.stderr)
        return None


def print_stats_table(stats_dict):
    """Formats and prints the statistics dictionary as a table."""
    if not stats_dict:
        print("No statistics to display.")
        return

    # Define column headers matching the user's image/request
    categories = ['total', 'iur', 'igr', '(?)', 'plus', 'minus', 'single', 'multi']
    # Format strings for header and rows {<LeftAlign:Width} {:>RightAlign:Width}
    header_format = "{:<12} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}"
    row_format = "{:<12} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}"

    print("\n" + "="*75)
    print("--- Summary Statistics per TFBS Family ---")
    print("="*75)
    print(header_format.format("TFBS_Family", *categories))
    print("-" * 75) # Adjust separator length based on format string widths

    # Sort families alphabetically for consistent output
    sorted_families = sorted(stats_dict.keys())

    for family in sorted_families:
        counts = stats_dict[family]
        # Use .get(key, 0) to safely access counts, defaulting to 0 if a key is missing
        print(row_format.format(
            family,
            counts.get('total', 0),
            counts.get('iur', 0),
            counts.get('igr', 0),
            counts.get('(?)', 0),
            counts.get('plus', 0),
            counts.get('minus', 0),
            counts.get('single', 0),
            counts.get('multi', 0)
        ))
    print("-" * 75 + "\n")


# --- Main Execution Block ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generates summary statistics from a multifasta file with tagged headers. "
                    "Assumes header format: '>coords tag_block number TFBS_Family [opt]'. "
                    "'single'/'multi' columns count based on the 'number' (1 vs 2-6).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input_fasta",
                        help="Path to the input tagged multifasta file.")
    # Optional output file for stats? For now, just print to console. Add later if needed.
    # parser.add_argument("-o", "--output", help="Optional path to save statistics table.")

    args = parser.parse_args()

    # Generate statistics by processing the input file
    statistics = generate_statistics(args.input_fasta)

    # Print the results table if statistics were successfully generated
    if statistics:
        print_stats_table(statistics)
    else:
        print("Statistics generation failed.", file=sys.stderr)
        sys.exit(1) # Exit with error status

    print("Script finished.", file=sys.stderr)
