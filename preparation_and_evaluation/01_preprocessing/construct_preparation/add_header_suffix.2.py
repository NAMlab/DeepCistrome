import argparse
import re

def rename_fasta_headers(input_list_file, multifasta_file, output_fasta_file):
    """
    Renames FASTA headers based on an input list file.

    Args:
        input_list_file (str): Path to the input list file.
        multifasta_file (str): Path to the input multifasta file.
        output_fasta_file (str): Path to the output fasta file with renamed headers.
    """

    header_map = {}
    try:
        with open(input_list_file, 'r') as infile_list:
            next(infile_list)  # Skip the header line
            for line in infile_list:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    identifier_range = parts[2]
                    new_header_suffix = f"{parts[0]}\t{parts[1]}"
                    header_map[identifier_range] = new_header_suffix
    except FileNotFoundError:
        print(f"Error: Input list file not found at {input_list_file}")
        return

    try:
        with open(multifasta_file, 'r') as infile_fasta, open(output_fasta_file, 'w') as outfile_fasta:
            header = None
            for line in infile_fasta:
                line = line.strip()
                if line.startswith('>'):
                    header = line
                    original_identifier = header.split(' ')[0][1:]  # Extract identifier after '>'

                    found_match = False
                    for identifier_range, suffix in header_map.items():
                        # Use regex to match the range at the beginning of the identifier
                        if re.match(r'^' + re.escape(identifier_range), original_identifier):
                            new_header = f"{header} {suffix}"
                            outfile_fasta.write(f"{new_header}\n")
                            found_match = True
                            break
                    if not found_match:
                        outfile_fasta.write(f"{header}\n")  # Write original header if no match
                else:
                    if header:
                        outfile_fasta.write(f"{line}\n")
    except FileNotFoundError:
        print(f"Error: Multifasta file not found at {multifasta_file}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rename FASTA headers based on an input list.")
    parser.add_argument("input_list", help="Path to the input list file.")
    parser.add_argument("multifasta", help="Path to the input multifasta file.")
    parser.add_argument("output_fasta", help="Path to the output fasta file.")

    args = parser.parse_args()

    rename_fasta_headers(args.input_list, args.multifasta, args.output_fasta)
    print(f"Renamed FASTA saved to: {args.output_fasta}")
