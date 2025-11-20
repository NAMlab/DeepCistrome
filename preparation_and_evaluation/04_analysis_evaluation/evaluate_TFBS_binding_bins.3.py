#!/usr/bin/env python3

import pandas as pd
import numpy as np
from sklearn.metrics import (
    recall_score,
    f1_score,
    matthews_corrcoef,
    average_precision_score,
)
from sklearn.utils import resample
from tqdm import tqdm
import warnings
import matplotlib.pyplot as plt
import os

# --- Configuration ---
# Define the base folder for inputs and outputs
EVALUATION_FOLDER = '../studies/deepCIS/deepCistome/EVALUATION/'

# Input and Output file paths
ACTUAL_FILE_PATH = os.path.join(EVALUATION_FOLDER, 'all_actual.csv')
PREDICTIONS_FILE_PATH = os.path.join(EVALUATION_FOLDER, 'all_predictions.csv')
OUTPUT_PLOT_FILE = os.path.join(EVALUATION_FOLDER, 'evaluation_performance_plot.png')
OUTPUT_TABLE_FILE = os.path.join(EVALUATION_FOLDER, 'evaluation_summary_table.csv')


# Number of bootstrap iterations for stable metrics
N_BOOTSTRAPS = 10

# Suppress warnings from sklearn for cases with no true labels in a bootstrap sample
warnings.filterwarnings("ignore", category=UserWarning)

def calculate_auprc_macro(y_true, y_pred):
    """
    Calculates macro-averaged AUPRC manually to avoid errors when a bootstrap
    sample for a specific class contains no positive instances.
    """
    auprc_per_class = []
    for i in range(y_true.shape[1]):
        y_true_class = y_true.iloc[:, i]
        if y_true_class.sum() > 0:
            y_pred_class = y_pred.iloc[:, i]
            auprc_per_class.append(average_precision_score(y_true_class, y_pred_class))
    return np.mean(auprc_per_class) if auprc_per_class else 0.0

def create_and_save_plot(results_df, n_bootstraps):
    """
    Generates and saves a joint bar plot with the counts plot first.
    """
    if results_df.empty:
        print("📊 Skipping plot generation as there are no results.")
        return

    print(f"\n📊 Generating final performance plot...")

    # --- Plotting Configuration ---
    metrics_to_plot = {
        'SENSITIVITY': 'Sensitivity (Micro Recall)',
        'F1_MICRO': 'F1 Score (Micro)',
        'MCC': 'Matthews Correlation Coefficient',
        'AUPRC_MICRO': 'AUPRC (Micro)'
    }
    grey_colors = ['#222222', '#555555', '#888888', '#BBBBBB']

    plot_data = results_df.copy()
    for metric in metrics_to_plot.keys():
        parts = plot_data[metric].str.split(' ± ', expand=True)
        plot_data[f'{metric}_mean'] = pd.to_numeric(parts[0])
        plot_data[f'{metric}_std'] = pd.to_numeric(parts[1])

    # Use a 3x2 grid to fit 1 counts plot + 4 performance plots
    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(8.27, 9.0))
    axes = axes.flatten()
    
    bins = plot_data['Bin (Predicted Sites)']

    # --- Generate Actual Counts Plot (1st plot) ---
    ax = axes[0]
    actual_counts = plot_data['Actual Counts']
    ax.bar(bins, actual_counts, color='#2c7fb8', alpha=0.8)

    ax.set_title("Total True Labels per Bin", fontsize=12)
    ax.set_ylabel("Sum of True Labels", fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', linestyle='--', alpha=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(axis='both', which='major', labelsize=9)

    # --- Generate Performance Metric Plots (4 plots, starting from the 2nd panel) ---
    for i, (metric_key, title) in enumerate(metrics_to_plot.items()):
        ax = axes[i + 1] # Start plotting on the second axis
        means = plot_data[f'{metric_key}_mean']
        stds = plot_data[f'{metric_key}_std']
        
        ax.bar(bins, means, yerr=stds, capsize=4, alpha=0.9, color=grey_colors[i])
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', linestyle='--', alpha=0.6)
        ax.set_axisbelow(True)

        ax.set_title(title, fontsize=12)
        ax.set_ylabel("Score", fontsize=10)
        ax.tick_params(axis='both', which='major', labelsize=9)
        
        y_min = 0
        y_max = max(1.0, (means + stds).max() * 1.05)
        ax.set_ylim(bottom=y_min, top=y_max)
    
    # --- Final Touches ---
    # Turn off the last, unused subplot
    axes[5].axis('off')

    # Set common X-axis labels for the bottom row panels
    for i in [2, 3]:
        axes[i].set_xlabel("Bin (Predicted Sites)", fontsize=11)
        if not bins.empty:
            tick_locations = np.arange(bins.min(), bins.max() + 1)
            axes[i].set_xticks(tick_locations)
            axes[i].set_xticklabels(tick_locations.astype(int))

    fig.suptitle(f'Model Performance by Number of Predicted Sites\n(n = {n_bootstraps} bootstraps)', fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    # Save the plot to the specified output folder
    plt.savefig(OUTPUT_PLOT_FILE, dpi=300, bbox_inches='tight')
    print(f"\n✅ Plot saved successfully to '{OUTPUT_PLOT_FILE}'")

def main():
    """
    Main function to load data, perform binned evaluation with bootstrapping,
    and save the results as a table and a plot to the EVALUATION folder.
    """
    print("🚀 Starting evaluation script...")
    
    os.makedirs(EVALUATION_FOLDER, exist_ok=True)
    
    try:
        y_true = pd.read_csv(ACTUAL_FILE_PATH, sep='\t')
        y_pred = pd.read_csv(PREDICTIONS_FILE_PATH, sep='\t')
    except FileNotFoundError as e:
        print(f"❌ Error: {e}. Please ensure input files exist in {EVALUATION_FOLDER}")
        return

    if y_true.shape != y_pred.shape:
        print("❌ Error: Shape mismatch between actual and prediction files.")
        return
    print(f"Data loaded successfully. Shape: {y_true.shape}")

    predicted_sites_count = y_pred.sum(axis=1)
    bins = sorted(predicted_sites_count.unique())
    print(f"Found {len(bins)} bins based on predicted site counts: {bins}")

    results_list = []
    for bin_val in bins:
        indices = predicted_sites_count[predicted_sites_count == bin_val].index
        y_true_bin = y_true.loc[indices]
        y_pred_bin = y_pred.loc[indices]
        n_samples = len(y_true_bin)

        if n_samples < 2:
            print(f"\nSkipping bin '{bin_val}' due to insufficient samples ({n_samples}).")
            continue

        print(f"\nProcessing Bin '{bin_val}' ({n_samples} samples)...")
        metric_scores = {
            'sensitivity': [], 'f1_macro': [], 'f1_micro': [],
            'mcc': [], 'auprc_macro': [], 'auprc_micro': []
        }

        for _ in tqdm(range(N_BOOTSTRAPS), desc=f"  Bootstrap Bin {bin_val}"):
            y_true_boot, y_pred_boot = resample(y_true_bin, y_pred_bin)
            if y_true_boot.empty: continue

            metric_scores['sensitivity'].append(recall_score(y_true_boot, y_pred_boot, average='micro', zero_division=0))
            metric_scores['f1_macro'].append(f1_score(y_true_boot, y_pred_boot, average='macro', zero_division=0))
            metric_scores['f1_micro'].append(f1_score(y_true_boot, y_pred_boot, average='micro', zero_division=0))
            metric_scores['mcc'].append(matthews_corrcoef(y_true_boot.values.ravel(), y_pred_boot.values.ravel()))
            metric_scores['auprc_micro'].append(average_precision_score(y_true_boot, y_pred_boot, average='micro'))
            metric_scores['auprc_macro'].append(calculate_auprc_macro(y_true_boot, y_pred_boot))

        bin_summary = {'Bin (Predicted Sites)': bin_val, 'Num Samples': n_samples}
        for name, scores in metric_scores.items():
            bin_summary[name.upper()] = f"{np.mean(scores):.3f} ± {np.std(scores):.3f}"
        
        bin_summary['Actual Counts'] = y_true_bin.values.sum()
        results_list.append(bin_summary)

    results_df = pd.DataFrame(results_list)
    print("\n\n--- 📋 Evaluation Table ---")
    if results_df.empty:
        print("No results were generated. Check bin sizes and data.")
    else:
        cols_order = [
            'Bin (Predicted Sites)', 'Num Samples', 'Actual Counts', 'SENSITIVITY', 
            'F1_MICRO', 'F1_MACRO', 'MCC', 'AUPRC_MICRO', 'AUPRC_MACRO'
        ]
        results_df = results_df[[col for col in cols_order if col in results_df.columns]]
        
        print(results_df.to_string())
        results_df.to_csv(OUTPUT_TABLE_FILE, index=False)
        print(f"\n✅ Summary table saved successfully to '{OUTPUT_TABLE_FILE}'")

    create_and_save_plot(results_df, N_BOOTSTRAPS)
    
    print("\n✅ Script finished.")

if __name__ == '__main__':
    main()
