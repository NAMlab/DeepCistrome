import sys
import numpy as np
from Bio import motifs, SeqIO

# --- CONFIGURATION ---
GENOME_FILE = "PATH/TO/GENOME.FA"  # Path to your Arabidopsis genome
# Replace with the path to your file containing the text you pasted
MOTIF_FILE = "JASPAR2026_CORE_plants_non-redundant_v2.jaspar" 
OUTPUT_FILE = "OUTPUT_TRACK.tsv"

# Threshold: Report hits with score > 85% of the max possible score
RELATIVE_SCORE_THRESHOLD = 0.85

def parse_jaspar_pfm_manually(filepath):

    print(f"[INFO] Parsing motifs from {filepath}...")
    
    parsed_motifs = []
    current_name = None
    counts = {}
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            
            # 1. Header Line (>MA0004.1 Arnt)
            if line.startswith(">"):
                # If we have a previous motif ready, save it
                if current_name and len(counts) == 4:
                    m = create_motif_from_counts(current_name, counts)
                    parsed_motifs.append(m)
                
                # Start new motif
                parts = line.split()
                # Use the ID (MA0004.1) as the name
                current_name = parts[0][1:] 
                counts = {} # Reset counts
            
            # 2. Data Lines (A [ 1 2 3 ])
            else:
                # Check which nucleotide this row belongs to
                base = line[0].upper()
                if base in ['A', 'C', 'G', 'T']:
                    # Remove brackets '[' and ']' and split by whitespace
                    content = line[1:].replace("[", "").replace("]", "")
                    # Convert number strings to integers
                    row_data = [float(x) for x in content.split()]
                    counts[base] = row_data

        # Save the very last motif
        if current_name and len(counts) == 4:
            m = create_motif_from_counts(current_name, counts)
            parsed_motifs.append(m)
            
    return parsed_motifs

def create_motif_from_counts(name, counts_dict):
    """Creates a Biopython Motif object from the raw dictionary."""
    # Ensure all rows have the same length
    length = len(counts_dict['A'])
    
    # Biopython expects a specific dictionary format
    m = motifs.Motif(alphabet="ACGT", counts=counts_dict)
    m.name = name
    return m

def main():
    # 1. Load Motifs
    try:
        motif_list = parse_jaspar_pfm_manually(MOTIF_FILE)
    except FileNotFoundError:
        print(f"[ERROR] Motif file not found: {MOTIF_FILE}")
        sys.exit(1)

    if not motif_list:
        print("[ERROR] No motifs parsed. Check input file.")
        sys.exit(1)

    print(f"[INFO] Loaded {len(motif_list)} motifs.")

    # 2. Pre-calculate PSSMs
    print("[INFO] Converting counts to Scoring Matrices (PSSMs)...")
    pssm_data = []
    
    for m in motif_list:
        # Normalize (Counts -> Probabilities) with pseudocounts
        pwm = m.counts.normalize(pseudocounts=0.5)
        
        # Convert to Log-Odds (PSSM)
        pssm = pwm.log_odds()
        
        # Calculate Threshold (85% of max possible score)
        max_score = pssm.max
        threshold = max_score * RELATIVE_SCORE_THRESHOLD
        
        pssm_data.append((m.name, pssm, threshold, len(m)))

    # 3. Scan Genome (Vectorized NumPy)
    print(f"[INFO] Scanning genome: {GENOME_FILE}")
    print(f"[INFO] Writing results to: {OUTPUT_FILE}")
    
    with open(OUTPUT_FILE, 'w') as out_f:
        for record in SeqIO.parse(GENOME_FILE, "fasta"):
            chrom_id = record.id
            seq_len = len(record.seq)
            print(f"  > Processing {chrom_id} ({seq_len} bp)...")
            
            # Biopython's pssm.calculate uses NumPy for speed
            seq = record.seq
            
            for name, pssm, thresh, m_len in pssm_data:
                
                # --- Forward Strand ---
                scores_fwd = pssm.calculate(seq)
                hits_fwd = np.where(scores_fwd >= thresh)[0]
                
                for pos in hits_fwd:
                    score = scores_fwd[pos]
                    out_f.write(f"{chrom_id}\t{pos}\t{pos + m_len}\t{name}\t{score:.2f}\t+\n")

                # --- Reverse Strand ---
                # Calculate on reverse complement
                rc_seq = seq.reverse_complement()
                scores_rev = pssm.calculate(rc_seq)
                hits_rev = np.where(scores_rev >= thresh)[0]
                
                for pos in hits_rev:
                    score = scores_rev[pos]
                    # Calculate genomic start coordinate
                    # Logic: SequenceLength - pos_on_rc - motif_length
                    gen_start = seq_len - pos - m_len
                    out_f.write(f"{chrom_id}\t{gen_start}\t{gen_start + m_len}\t{name}\t{score:.2f}\t-\n")

    print("[SUCCESS] Done.")

if __name__ == "__main__":
    main()
