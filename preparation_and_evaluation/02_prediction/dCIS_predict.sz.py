# dCIS_predict_hardcoded_tflist.py

import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
import os
import sys
from tensorflow.keras.models import load_model
from typing import Any # Needed for embedded one_hot_encode type hints

# Make sure utils.py containing one_hot_encode is in the same directory
# or accessible in the Python path.
try:
    from utils import one_hot_encode
except ImportError:
    print("Info: 'utils.py' not found. Using embedded 'one_hot_encode' function.", file=sys.stderr)
    # Define it here if utils.py is unavailable, using the code you provided:
    def one_hot_encode(sequence: str,
                       alphabet: str = 'ACGT',
                       neutral_alphabet: str = 'N',
                       neutral_value: Any = 0,
                       dtype=np.float32) -> np.ndarray:
        """One-hot encode sequence."""
        def to_uint8(string):
            return np.frombuffer(string.encode('ascii'), dtype=np.uint8)
        hash_table = np.zeros((np.iinfo(np.uint8).max, len(alphabet)), dtype=dtype)
        hash_table[to_uint8(alphabet)] = np.eye(len(alphabet), dtype=dtype)
        hash_table[to_uint8(neutral_alphabet)] = neutral_value
        hash_table = hash_table.astype(dtype)
        # Handle potential non-ACGTN characters gracefully - map them to neutral value
        # Create a buffer for the sequence to handle potential errors
        seq_bytes = to_uint8(sequence.upper()) # Convert to uppercase for consistency
        # Ensure all bytes are within the hash_table range, map unknowns to 'N' ASCII value (78)
        valid_chars = to_uint8(alphabet + neutral_alphabet)
        # Create mask of invalid characters
        invalid_mask = ~np.isin(seq_bytes, valid_chars)
        # Map invalid characters to 'N' (ASCII 78)
        seq_bytes[invalid_mask] = 78
        return hash_table[seq_bytes]

# --- Configuration ---
# Removed TF_FAM_INFO_PATH - using hardcoded list now
MODEL_PATH_DEFAULT = 'saved_models/model_chrom_1_model.h5'
EXPECTED_LENGTH_DEFAULT = 250 # Model expects sequences of this length
PREDICTION_THRESHOLD_DEFAULT = 0.5
OUTPUT_SEPARATOR = '\t' # Use tab for the output matrix

# Hardcoded list of 46 TF names corresponding to model output channels
TF_NAMES_RAW = "ABI3VP1_tnt AP2EREBP_tnt ARF_ecoli ARF_tnt ARID_tnt BBRBPC_tnt BES1_tnt BZR_tnt C2C2YABBY_tnt C2C2dof_tnt C2C2gata_tnt C2H2_tnt C3H_tnt CAMTA_tnt CPP_tnt E2FDP_tnt EIL_tnt FAR1_tnt G2like_tnt GRF_tnt GeBP_tnt HB_tnt HSF_tnt Homeobox_ecoli Homeobox_tnt LOBAS2_tnt MADS_tnt MYB_tnt MYBrelated_tnt NAC_tnt ND_tnt Orphan_tnt RAV_tnt REM_tnt RWPRK_tnt S1Falike_tnt SBP_tnt SRS_tnt TCP_tnt Trihelix_tnt WRKY_tnt ZFHD_tnt bHLH_tnt bZIP_tnt mTERF_tnt zfGRF_tnt"

# Process the hardcoded list (split and clean)
tf_names_list = TF_NAMES_RAW.split()
# Apply the same cleaning step as before for consistency, if desired
tf_names_list = [name.replace('_tnt', '') for name in tf_names_list]

# --- End Configuration ---

pd.options.display.width = 0

# --- Argument Parsing ---
parser = argparse.ArgumentParser(
    description="Predict TF binding using a Keras model on sequences from a FASTA file."
)
parser.add_argument("input_fasta", help="Path to the input FASTA file.")
parser.add_argument(
    "-o", "--output_csv",
    required=True,
    help="Path for the output predictions TSV (tab-separated values) file."
)
parser.add_argument(
    "--model",
    default=MODEL_PATH_DEFAULT,
    help=f"Path to the trained Keras model file (default: {MODEL_PATH_DEFAULT})"
)
# Removed --tf_list argument
parser.add_argument(
    "--length",
    type=int,
    default=EXPECTED_LENGTH_DEFAULT,
    help=f"Expected sequence length for the model (default: {EXPECTED_LENGTH_DEFAULT})"
)
parser.add_argument(
    "--threshold",
    type=float,
    default=PREDICTION_THRESHOLD_DEFAULT,
    help=f"Prediction threshold (default: {PREDICTION_THRESHOLD_DEFAULT})"
)

args = parser.parse_args()
# Update config based on args
MODEL_PATH = args.model
EXPECTED_LENGTH = args.length
PREDICTION_THRESHOLD = args.threshold
# --- End Argument Parsing ---

# Verify the hardcoded list has 46 names
N_EXPECTED_TFS = 46
if len(tf_names_list) != N_EXPECTED_TFS:
    print(f"Error: The hardcoded TF list contains {len(tf_names_list)} names, but {N_EXPECTED_TFS} were expected.", file=sys.stderr)
    print(f"Hardcoded list used: {tf_names_list}", file=sys.stderr)
    sys.exit(1)
else:
    print(f"Using hardcoded list of {len(tf_names_list)} TF names.")
    # print(f" TF names: {tf_names_list}") # Optional: print the list

# Load the Keras model
try:
    print(f"Loading model from {MODEL_PATH}...")
    model = load_model(MODEL_PATH)
    print("Model loaded successfully.")
    # Optional: Check model input/output shape if possible
    try:
        model_input_shape = model.input_shape
        model_output_shape = model.output_shape
        print(f"  Model expected input shape: {model_input_shape}, output shape: {model_output_shape}")
        if len(model_input_shape) != 3 or model_input_shape[1] != EXPECTED_LENGTH or model_input_shape[2] != 4:
             print(f"  Warning: Model input shape {model_input_shape} might not match expected length {EXPECTED_LENGTH} and 4 channels.", file=sys.stderr)
        if len(model_output_shape) != 2 or model_output_shape[1] != N_EXPECTED_TFS:
             print(f"  Warning: Model output shape {model_output_shape} might not match expected number of TFs ({N_EXPECTED_TFS}).", file=sys.stderr)
    except Exception as e:
        print(f"  Info: Could not verify model input/output shape automatically ({e}).", file=sys.stderr)

except FileNotFoundError:
    print(f"Error: Model file not found at '{MODEL_PATH}'", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"Error loading Keras model: {e}", file=sys.stderr)
    sys.exit(1)


# --- Process specified FASTA file ---
input_filepath = args.input_fasta
print(f"\n--- Processing file: {input_filepath} ---")

seqs_enc, seq_full_headers = [], []
seq_count = 0
valid_seq_count = 0

try:
    # Read sequences and filter by length
    for rec in SeqIO.parse(handle=input_filepath, format='fasta'):
        seq_count += 1
        # Check if sequence has the exact expected length
        if len(rec.seq) == EXPECTED_LENGTH:
            valid_seq_count += 1
            try:
                # One-hot encode the sequence
                encoded = one_hot_encode(str(rec.seq))
                seqs_enc.append(encoded)
                # Store the full header line, using rec.description
                seq_full_headers.append(rec.description)
            except Exception as e:
                print(f"  Warning: Error encoding sequence {rec.id}: {e}", file=sys.stderr)
        # else: # Optional: uncomment to log sequences with wrong length
            # print(f"  Skipping sequence {rec.id}: Length {len(rec.seq)} != {EXPECTED_LENGTH}", file=sys.stderr)

except FileNotFoundError:
    print(f"Error: Input FASTA file not found at '{input_filepath}'", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"Error reading FASTA file '{input_filepath}': {e}", file=sys.stderr)
    sys.exit(1)


print(f"  Found {seq_count} sequences, {valid_seq_count} have the expected length ({EXPECTED_LENGTH} bp).")

# Proceed only if there are valid sequences to predict on
if not seqs_enc:
    print("Error: No sequences of the expected length found in the input file. Cannot generate predictions.", file=sys.stderr)
    sys.exit(1)

# Convert lists to NumPy arrays
seqs_enc = np.array(seqs_enc)
# seq_full_headers remains a list of strings

print(f"  Encoded sequences shape: {seqs_enc.shape}, Number of headers: {len(seq_full_headers)}")

# Perform prediction directly
try:
    print("  Running model predictions...")
    predictions_raw = model.predict(seqs_enc)

    # Apply threshold
    predictions_bool = predictions_raw > PREDICTION_THRESHOLD
    predictions_int = predictions_bool.astype('int')

    # Shape should now be (num_valid_seqs, num_tf_families)
    print(f"  Predictions shape after thresholding: {predictions_int.shape}")

    # Verify prediction output size against the hardcoded list size
    if predictions_int.shape[0] != len(seq_full_headers):
        print(f"  Error: Number of predictions ({predictions_int.shape[0]}) does not match number of valid sequences ({len(seq_full_headers)}). Internal error.", file=sys.stderr)
        sys.exit(1)

    if predictions_int.shape[1] != len(tf_names_list):
         print(f"\n--- ERROR ---", file=sys.stderr)
         print(f"Mismatch between model output size and hardcoded TF list size:", file=sys.stderr)
         print(f"  Model produced {predictions_int.shape[1]} outputs per sequence.", file=sys.stderr)
         print(f"  The hardcoded list contains {len(tf_names_list)} TF names.", file=sys.stderr)
         print(f"  Please check the hardcoded 'TF_NAMES_RAW' list in the script or the model '{MODEL_PATH}'.", file=sys.stderr)
         print(f"--- END ERROR ---", file=sys.stderr)
         sys.exit(1) # Exit script cleanly

    # Create results DataFrame using the hardcoded TF names
    predictions_df = pd.DataFrame(predictions_int, columns=tf_names_list)

    # Insert the full headers as the first column
    predictions_df.insert(0, 'SequenceHeader', seq_full_headers)

    print(f"  Generated predictions DataFrame shape: {predictions_df.shape}")
    print("  Predictions head:\n", predictions_df.head())

    # Save predictions to CSV (using specified separator, no index)
    output_filepath = args.output_csv
    predictions_df.to_csv(output_filepath, sep=OUTPUT_SEPARATOR, index=False)
    print(f"  Predictions saved to: {output_filepath}")

except Exception as e:
    print(f"  Error during prediction or saving results: {e}", file=sys.stderr)
    sys.exit(1)

print("\nScript finished successfully.")
