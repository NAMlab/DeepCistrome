#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os

# --- Configuration ---
# Define the base folder where the input file is and outputs will be saved
EVALUATION_FOLDER = '../studies/deepCIS/deepCistome/EVALUATION/'

# Input and Output file paths
PREDICTIONS_FILE_PATH = os.path.join(EVALUATION_FOLDER, 'all_predictions.csv')
OUTPUT_FILE_PATH = os.path.join(EVALUATION_FOLDER, 'average_cooccurrence_bin.csv')


def calculate_average_bin(predictions_df):
    """
    For each TF, calculates the average number of total co-occurring TFs.

    Args:
        predictions_df (pd.DataFrame): DataFrame with binary predictions.

    Returns:
        pd.Series: A Series with TF families as the index and their
                   average co-occurrence bin as the value.
    """
    print("Calculating the total number of predicted sites for each sample...")
    # Calculate the bin number for each sample (row sum)
    bin_assignments = predictions_df.sum(axis=1)

    # Dictionary to store the single value for each TF
    avg_bin_results = {}

    print("Calculating average co-occurrence bin for each TF family...")
    # Iterate over each TF family (column)
    for tf_family in predictions_df.columns:
        # Find all samples where this TF is present (value is 1)
        tf_is_present = predictions_df[tf_family] == 1

        if tf_is_present.sum() == 0:
            # If a TF is never predicted, its average bin is undefined (or 0)
            avg_bin = 0.0
        else:
            # For the samples where the TF is present, get their bin numbers
            bins_for_this_tf = bin_assignments[tf_is_present]
            
            # Calculate the mean of these bin numbers
            avg_bin = bins_for_this_tf.mean()
        
        avg_bin_results[tf_family] = avg_bin

    # Convert dictionary to a pandas Series for easy sorting and display
    results_series = pd.Series(avg_bin_results)
    
    # Sort the results to easily see which TFs have the highest co-occurrence
    results_series = results_series.sort_values(ascending=False)
    
    return results_series


def main():
    """
    Main function to load data, run calculation, and save results.
    """
    print("🚀 Starting average co-occurrence bin calculation...")
    
    os.makedirs(EVALUATION_FOLDER, exist_ok=True)
    
    try:
        print(f"Loading predictions from: {PREDICTIONS_FILE_PATH}")
        predictions_df = pd.read_csv(PREDICTIONS_FILE_PATH, sep='\t')
    except FileNotFoundError as e:
        print(f"❌ Error: {e}. Please ensure the input file exists.")
        return
        
    print(f"Data loaded successfully. Shape: {predictions_df.shape}")

    # Perform the calculation
    avg_bin_series = calculate_average_bin(predictions_df)
    
    # Convert Series to DataFrame for saving with a proper header
    results_df = avg_bin_series.to_frame(name='Average_Cooccurrence_Bin')
    results_df.index.name = 'TF_Family'

    try:
        results_df.to_csv(OUTPUT_FILE_PATH)
        print(f"\n✅ Results saved successfully to '{OUTPUT_FILE_PATH}'")
    except Exception as e:
        print(f"\n❌ Error saving file: {e}")

    # Display the results table
    print("\n--- Average Co-occurrence Bin per TF Family ---")
    print(results_df)
    print("\n✅ Script finished.")


if __name__ == '__main__':
    main()
