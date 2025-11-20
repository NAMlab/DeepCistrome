#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import sys
import os
import argparse

# --- Function Definitions ---

def parse_gff(gff_filepath):
    """
    Parses a GFF3 file to extract gene information (chr, start, end, strand).

    Args:
        gff_filepath (str): Path to the GFF3 file.

    Returns:
        dict: A dictionary where keys are chromosome names and values are lists
              of tuples, each tuple containing (gene_start, gene_end, strand).
              Returns None if the file cannot be opened or an error occurs.
    """
    genes_by_chr = {}
    try:
        with open(gff_filepath, 'r') as gff_file:
            for line_num, line in enumerate(gff_file, 1):
                if line.startswith('#'):
                    continue # Skip comment lines
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue # Skip lines that don't have enough fields

                feature_type = fields[2]
                # Focus only on 'gene' features (case-insensitive check)
                if feature_type.lower() == 'gene':
                    chromosome = fields[0]
                    try:
                        # GFF is 1-based, inclusive start/end
                        start = int(fields[3])
                        end = int(fields[4])
                    except ValueError:
                        print(f"Warning: Skipping GFF line {line_num} due to invalid coordinates: {line.strip()}", file=sys.stderr)
                        continue

                    strand = fields[6]
                    if strand not in ['+', '-']:
                        print(f"Warning: Skipping GFF line {line_num} with unknown strand '{strand}': {line.strip()}", file=sys.stderr)
                        continue

                    # Store gene info
                    if chromosome not in genes_by_chr:
                        genes_by_chr[chromosome] = []
                    genes_by_chr[chromosome].append((start, end, strand))

    except IOError as e:
        print(f"Error reading GFF file {gff_filepath}: {e}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred during GFF parsing: {e}", file=sys.stderr)
        return None

    print(f"Finished parsing GFF. Found genes on {len(genes_by_chr)} chromosomes.", file=sys.stderr)
    # Optional: Sort genes by start position for potential minor optimization later
    for chrom in genes_by_chr:
        genes_by_chr[chrom].sort()
    return genes_by_chr


def process_entry(header, sequence, genes_by_chr, iur_distance, header_pattern):
    """
    Processes a single fasta entry: finds iur/igr region tag and orientation
    based on whether the sequence window's CENTER POINT falls within the
    defined regions (iur: upstream promoter, igr: within gene body).
    iur tags take priority over igr tags.
    """
    match = header_pattern.match(header)
    final_tag = "(?)" # Default tag if no match found or coords invalid

    if match:
        seq_chr = match.group(1)
        center_point = None # Initialize center point
        seq_coords_valid = False # Flag for coordinate validity

        try:
            # Fasta coords often derived from GFF, assume 1-based inclusive
            seq_start = int(match.group(2))
            seq_end = int(match.group(3))
            if seq_start > seq_end:
                 print(f"Warning: Sequence start ({seq_start}) > end ({seq_end}) in header: {header}", file=sys.stderr)
                 # Coordinates are invalid, cannot calculate center point
            else:
                 # Calculate center point (using float division)
                 center_point = (seq_start + seq_end) / 2.0
                 seq_coords_valid = True
        except ValueError:
             print(f"Warning: Could not parse coordinates from header: {header}", file=sys.stderr)
             # Coordinates are invalid

        rest_of_header = match.group(4).strip()

        # Proceed only if coordinates were valid, center_point calculated, and chromosome exists in GFF
        if seq_coords_valid and center_point is not None and seq_chr in genes_by_chr:

            # --- First Pass: Check if CENTER POINT falls in IUR (Promoter) Region ---
            for gene_start, gene_end, gene_strand in genes_by_chr[seq_chr]:
                is_iur_match = False
                if gene_strand == '+':
                    # IUR region: [gene_start - iur_distance, gene_start - 1]
                    promoter_start = gene_start - iur_distance
                    promoter_end = gene_start - 1
                    # Check if center point is within the region
                    if center_point >= promoter_start and center_point <= promoter_end:
                        is_iur_match = True
                elif gene_strand == '-':
                    # IUR region: [gene_end + 1, gene_end + iur_distance]
                    promoter_start = gene_end + 1
                    promoter_end = gene_end + iur_distance
                    # Check if center point is within the region
                    if center_point >= promoter_start and center_point <= promoter_end:
                        is_iur_match = True

                if is_iur_match:
                    final_tag = f"iur({gene_strand})"
                    break # IUR found, highest priority, stop searching

            # --- Second Pass: Check if CENTER POINT falls in IGR (IntraGenic Region) (only if no IUR was found) ---
            if final_tag == "(?)": # Only proceed if no IUR match was found
                for gene_start, gene_end, gene_strand in genes_by_chr[seq_chr]:
                    is_igr_match = False
                    # Check if center point is within the gene boundaries [start, end]
                    # GFF gene start/end are inclusive
                    if center_point >= gene_start and center_point <= gene_end:
                         is_igr_match = True

                    if is_igr_match:
                        # Assign tag based on the strand of the gene it falls within
                        final_tag = f"igr({gene_strand})"
                        break # IGR found, stop searching

            # --- Construct the new header ---
            coords_part = f"{seq_chr}:{seq_start}-{seq_end}"
            if rest_of_header:
                 new_header = f">{coords_part} {final_tag} {rest_of_header}"
            else:
                 new_header = f">{coords_part} {final_tag}"

        # --- Handle cases where classification wasn't possible ---
        elif seq_coords_valid and seq_chr not in genes_by_chr:
             # Coords valid, but Chr not in GFF -> Use default "(?)" tag
             coords_part = f"{seq_chr}:{seq_start}-{seq_end}"
             print(f"Warning: Chromosome '{seq_chr}' from header '{header}' not found in GFF data.", file=sys.stderr)
             if rest_of_header:
                 new_header = f">{coords_part} {final_tag} {rest_of_header}" # final_tag is "(?)"
             else:
                 new_header = f">{coords_part} {final_tag}" # final_tag is "(?)"
        else:
            # Coordinates invalid/unparseable from header -> Use default "(?)" tag on original header
             new_header = f"{header} {final_tag}" # final_tag is "(?)"

    else: # Header format doesn't match expected pattern at all
        print(f"Warning: Header format not recognized, cannot parse coordinates: {header}", file=sys.stderr)
        new_header = f"{header} {final_tag}" # Append default marker "(?)"

    # --- Print Output ---
    # Print modified header
    print(new_header)
    # Print sequence, wrapped at 60 characters per line
    line_length = 60
    for i in range(0, len(sequence), line_length):
        print(sequence[i:i+line_length])


def process_fasta(fasta_filepath, genes_by_chr, iur_upstream_distance):
    """
    Reads a multifasta file, determines orientation and region type (iur/igr)
    for each entry based on gene proximity using its center point,
    and prints the modified fasta entries.
    """
    if genes_by_chr is None:
        print("Error: Gene data is not available. Cannot process fasta.", file=sys.stderr)
        return

    header = None
    sequence_parts = []
    # Regex to capture chromosome (non-whitespace), start (digits), end (digits)
    # Allows for optional extra info after a space
    header_pattern = re.compile(r"^>(\S+):(\d+)-(\d+)\s*(.*)")

    try:
        with open(fasta_filepath, 'r') as fasta_file:
            for line_num, line in enumerate(fasta_file, 1):
                line = line.strip()
                if not line:
                    continue # Skip empty lines

                if line.startswith('>'):
                    # Process the previous entry if one exists
                    if header:
                        # Call process_entry for the completed entry
                        process_entry(header, "".join(sequence_parts), genes_by_chr, iur_upstream_distance, header_pattern)
                        sequence_parts = [] # Reset sequence accumulator

                    header = line # Store the new header line
                else:
                    # Append sequence line if we are currently processing an entry
                    if header:
                       sequence_parts.append(line)
                    # Optional: could add sequence character validation here

            # Process the very last entry in the file after the loop ends
            if header:
                 process_entry(header, "".join(sequence_parts), genes_by_chr, iur_upstream_distance, header_pattern)

    except IOError as e:
        print(f"Error reading Fasta file {fasta_filepath}: {e}", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred during Fasta processing: {e}", file=sys.stderr)

    print(f"Finished processing Fasta.", file=sys.stderr)


# --- Main Execution ---
if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Adds orientation (+/-) and region type (iur/igr) tags to fasta headers "
                    "based on the location of the sequence window's center point relative to "
                    "genes defined in a GFF3 file. 'iur' (upstream promoter) tags take priority "
                    "over 'igr' (within gene body) tags. Unclassified regions are marked '(?)'.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter # Shows default values in help
    )

    # Required positional arguments
    parser.add_argument("gff_file",
                        help="Path to the input GFF3 annotation file.")
    parser.add_argument("fasta_file",
                        help="Path to the input multifasta sequence file.")

    # Optional arguments
    parser.add_argument("-d", "--distance", type=int, default=1000, dest="iur_distance",
                        help="Upstream distance (bp) from gene TSS to define the 'iur' (promoter) region.")
    parser.add_argument("-o", "--output",
                         help="Optional path to write the output fasta file. "
                              "If not specified, output goes to standard output (console).")


    # Parse arguments from command line
    args = parser.parse_args()

    # Assign arguments to variables
    gff_filename = args.gff_file
    fasta_filename = args.fasta_file
    promoter_distance_for_iur = args.iur_distance # Use the parsed value for IUR distance
    output_filename = args.output

    # --- Input File Validation ---
    if not os.path.isfile(gff_filename):
        print(f"Error: GFF file not found or is not a file: {gff_filename}", file=sys.stderr)
        sys.exit(1) # Exit with a non-zero code indicating an error
    if not os.path.isfile(fasta_filename):
        print(f"Error: Fasta file not found or is not a file: {fasta_filename}", file=sys.stderr)
        sys.exit(1)

    # --- Start Processing ---
    # Print status messages to standard error
    print(f"Using GFF file: {gff_filename}", file=sys.stderr)
    print(f"Using Fasta file: {fasta_filename}", file=sys.stderr)
    print(f"Using 'iur' promoter distance: {promoter_distance_for_iur} bp", file=sys.stderr)
    print(f"Checking for 'igr' (intragenic) regions based on center point.", file=sys.stderr)
    if output_filename:
        print(f"Output will be written to: {output_filename}", file=sys.stderr)
    else:
        print("Output will be written to standard output.", file=sys.stderr)

    # 1. Parse the GFF file to get gene data
    gene_data = parse_gff(gff_filename)

    # 2. Process the Fasta file using the gene data
    if gene_data is not None: # Proceed only if GFF parsing was successful
        original_stdout = sys.stdout # Keep track of the original standard output
        output_stream = None
        # Redirect stdout if an output file is specified
        if output_filename:
            try:
                output_stream = open(output_filename, 'w')
                sys.stdout = output_stream # Point standard output to the file
            except IOError as e:
                # If output file cannot be opened, print error to stderr and exit
                print(f"Error: Could not open output file {output_filename} for writing: {e}", file=sys.stderr)
                sys.exit(1)

        # Call the main processing function - output goes to sys.stdout (file or console)
        process_fasta(fasta_filename, gene_data, promoter_distance_for_iur)

        # Close the output file and restore standard output if it was redirected
        if output_stream:
            output_stream.close()
            sys.stdout = original_stdout # Restore standard output
            print(f"Output successfully written to {output_filename}", file=sys.stderr)

    else:
        # GFF parsing failed earlier
        print("\nExiting due to errors in GFF processing.", file=sys.stderr)
        sys.exit(1) # Exit with an error code

    print("\nScript finished.", file=sys.stderr)
# --------------------
