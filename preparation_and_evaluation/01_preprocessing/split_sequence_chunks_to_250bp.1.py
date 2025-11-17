#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
import sys
import textwrap # Although SeqIO handles wrapping, good practice to import if needed

def split_fasta_overlap(input_fasta, output_fasta, chunk_size=250, overlap=50):
    """
    Splits sequences in a FASTA file into smaller, potentially overlapping chunks.

    Args:
        input_fasta (str): Path to the input FASTA file.
        output_fasta (str): Path to the output FASTA file for chunked sequences.
        chunk_size (int): The desired length of each chunk.
        overlap (int): The number of base pairs overlap between consecutive chunks.
    """
    if overlap >= chunk_size:
        print(f"Error: Overlap ({overlap}) cannot be greater than or equal to chunk size ({chunk_size}).", file=sys.stderr)
        sys.exit(1)
    if overlap < 0:
        print(f"Error: Overlap ({overlap}) cannot be negative.", file=sys.stderr)
        sys.exit(1)

    step_size = chunk_size - overlap
    if step_size <= 0:
         print(f"Error: Step size (chunk_size - overlap = {step_size}) must be positive.", file=sys.stderr)
         sys.exit(1)


    print(f"Reading sequences from: {input_fasta}")
    print(f"Splitting into chunks of size: {chunk_size} bp")
    print(f"Overlapping by: {overlap} bp (step size: {step_size} bp)")

    output_records = []
    total_seqs_processed = 0
    total_chunks_generated = 0

    try:
        with open(input_fasta, 'r') as infile:
            for record in SeqIO.parse(infile, "fasta"):
                total_seqs_processed += 1
                original_seq = record.seq
                original_id = record.description # Use full description for context
                seq_len = len(original_seq)

                if seq_len < chunk_size:
                    print(f"  Warning: Sequence '{record.id}' (length {seq_len}) is shorter than chunk size ({chunk_size}). Skipping.", file=sys.stderr)
                    continue

                chunk_index = 1
                start_0based = 0
                # Loop as long as a *full* chunk can be extracted from the current start position
                while start_0based + chunk_size <= seq_len:
                    end_0based = start_0based + chunk_size
                    chunk_seq_str = str(original_seq[start_0based:end_0based])

                    # Create a new ID/description for the chunk
                    # Format: OriginalDescription_chunkN_start-end (1-based coords within original)
                    chunk_id = f"{original_id}_chunk{chunk_index}_{start_0based+1}-{end_0based}"
                    chunk_record = SeqRecord(Seq(chunk_seq_str), id=chunk_id, description="") # ID has all info

                    output_records.append(chunk_record)
                    total_chunks_generated += 1
                    chunk_index += 1

                    # Move the start position for the next chunk
                    start_0based += step_size

                # Note: This loop intentionally only creates full-sized chunks.
                # Any remaining sequence at the end shorter than chunk_size is ignored.

                if total_seqs_processed % 1000 == 0:
                     print(f"  Processed {total_seqs_processed} original sequences...")


    except FileNotFoundError:
        print(f"Error: Input FASTA file not found at '{input_fasta}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while reading '{input_fasta}': {e}", file=sys.stderr)
        sys.exit(1)

    if not output_records:
        print("Error: No chunks were generated. Check input file, chunk size, and overlap.", file=sys.stderr)
        sys.exit(1)

    print(f"\nProcessed {total_seqs_processed} original sequences.")
    print(f"Generated {total_chunks_generated} chunk sequences.")

    try:
        with open(output_fasta, 'w') as outfile:
            # Use SeqIO.write for proper FASTA formatting
            SeqIO.write(output_records, outfile, "fasta")
        print(f"Chunked sequences written to: {output_fasta}")
    except IOError as e:
        print(f"Error writing to output file '{output_fasta}': {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split sequences in a FASTA file into potentially overlapping chunks of a specified size.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input_fasta", help="Path to the input FASTA file.")
    parser.add_argument("output_fasta", help="Path for the output FASTA file with chunked sequences.")
    parser.add_argument("-c", "--chunk_size", type=int, default=250,
                        help="The desired length of each sequence chunk.")
    parser.add_argument("-ov", "--overlap", type=int, default=50, # Defaulting overlap to 50 as requested
                        help="Number of base pairs of overlap between consecutive chunks.")

    args = parser.parse_args()

    # Additional validation is now inside the main function
    split_fasta_overlap(args.input_fasta, args.output_fasta, args.chunk_size, args.overlap)
