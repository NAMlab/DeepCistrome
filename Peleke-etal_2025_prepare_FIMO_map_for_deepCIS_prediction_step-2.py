import sys
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# --- CONFIGURATION ---
GENOME_FILE = "PATH/TO/GENOME.FA"             # Your Arabidopsis Genome
BED_FILE = "AthaTAIR_dCISipm_tfbs_hits.bed"                     # Your results from the previous step
OUTPUT_FASTA = "OUTPUT_TRACK.bed"   # Output file for the prediction tool
TARGET_LEN = 250                               # Required length for DeepCIS

def load_genome(fasta_path):
    """Loads the genome into a dictionary for fast access."""
    print(f"[INFO] Loading genome from {fasta_path}...")
    genome = {}
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            genome[record.id] = record.seq
        print(f"[INFO] Loaded {len(genome)} chromosomes.")
        return genome
    except FileNotFoundError:
        print(f"[ERROR] Genome file not found: {fasta_path}")
        sys.exit(1)

def extract_and_pad(seq_obj, start, end, strand, chrom_len):
    """
    Extracts sequence. Handles boundaries by padding with N if needed.
    Returns (sequence_string, actual_start, actual_end)
    """
    # 0-based coordinates for slicing
    slice_start = max(0, start)
    slice_end = min(chrom_len, end)
    
    # Extract valid DNA
    seq_slice = seq_obj[slice_start:slice_end]
    
    # Pad Left (if we tried to extract before start of chromosome)
    pad_left = 0
    if start < 0:
        pad_left = abs(start)
    
    # Pad Right (if we tried to extract past end of chromosome)
    pad_right = 0
    if end > chrom_len:
        pad_right = end - chrom_len
        
    # Construct final sequence with N padding
    final_seq = ("N" * pad_left) + str(seq_slice) + ("N" * pad_right)
    
    # Handle Strand (Reverse Complement if needed)
    if strand == '-':
        final_seq = str(Seq(final_seq).reverse_complement())
        
    return final_seq

def main():
    # 1. Load Genome
    genome = load_genome(GENOME_FILE)
    
    # 2. Read BED file
    print(f"[INFO] Reading BED file: {BED_FILE}")
    try:
        # Assuming standard BED columns: chrom, start, end, name, score, strand
        cols = ['chrom', 'start', 'end', 'name', 'score', 'strand']
        df = pd.read_csv(BED_FILE, sep='\t', header=None, names=cols, usecols=[0,1,2,3,4,5])
    except Exception as e:
        print(f"[ERROR] Could not read BED file: {e}")
        sys.exit(1)

    print(f"[INFO] Processing {len(df)} motif hits...")
    
    with open(OUTPUT_FASTA, 'w') as out_f:
        count = 0
        
        for _, row in df.iterrows():
            chrom = str(row['chrom'])
            motif_start = int(row['start'])
            motif_end = int(row['end'])
            motif_id = row['name']
            strand = row['strand']
            
            if chrom not in genome:
                continue

            # --- Calculate Center and Window ---
            # Calculate integer center of the motif
            motif_center = (motif_start + motif_end) // 2
            
            # Determine window start/end (half of 250 is 125)
            # Since Python slicing is [start:end), we go -125 and +125
            win_start = motif_center - (TARGET_LEN // 2)
            win_end = win_start + TARGET_LEN
            
            # --- Extract Sequence ---
            # Retrieve chromosome sequence object
            chrom_seq = genome[chrom]
            chrom_len = len(chrom_seq)
            
            sequence = extract_and_pad(chrom_seq, win_start, win_end, strand, chrom_len)
            
            # Double check length
            if len(sequence) != TARGET_LEN:
                print(f"[WARNING] Generated sequence length {len(sequence)} != 250. Skipping.")
                continue

            # --- Generate ID ---
            # Format: MotifID_Chrom:MotifStart-MotifEnd_Window:Start-End(Strand)
            header_id = f"{motif_id}_{chrom}:{motif_start}-{motif_end}_win:{win_start}-{win_end}({strand})"
            
            # Write to FASTA
            out_f.write(f">{header_id}\n{sequence}\n")
            count += 1
            
            if count % 50000 == 0:
                print(f"  > Processed {count} sequences...")

    print(f"[SUCCESS] Extracted {count} sequences to {OUTPUT_FASTA}")
    print("You can now run the deepCIS prediction tool.")

if __name__ == "__main__":
    main()
