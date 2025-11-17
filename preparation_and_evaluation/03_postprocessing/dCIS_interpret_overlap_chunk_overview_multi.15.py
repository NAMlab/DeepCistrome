#!/usr/bin/env python3
import argparse
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# --- Constants ---
TOTAL_SEQUENCE_LENGTH = 1500
WINDOW_SIZE = 250
STEP_SIZE = 10
MAX_WINDOW_OVERLAP = WINDOW_SIZE // STEP_SIZE

# --- Plotting Enhancements for Nature Publication Standards ---
DPI = 600
TITLE_FONTSIZE = 20
LABEL_FONTSIZE = 16
TICK_FONTSIZE = 12
# Using a colorblind-friendly and professional color palette
NATURE_PALETTE = {
    'count': '#377eb8',      # Blue
    'pi': '#ff7f00',       # Orange
    'proportion': '#4daf4a',   # Green
    'uncertainty': '#a65628', # Brown (New)
    'start': '#e41a1c',       # Red
    'end': '#984ea3'        # Purple
}
PALETTE_START_END = {'Start': NATURE_PALETTE['start'], 'End': NATURE_PALETTE['end']}

def parse_identifier(identifier_str):
    """Parses an identifier string to extract the base sequence identifier and chunk number."""
    match = re.search(r'(.*)_chunk(\d+)', str(identifier_str))
    if match:
        return match.group(1), int(match.group(2))
    return None, None

def analyze_and_summarize_events(input_file, output_dir):
    """
    Main function to detect TFBS events, calculate summary statistics, and generate plots.
    """
    print(f"Loading data from: {input_file}")
    input_path = Path(input_file)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    try:
        df = pd.read_csv(input_path, header=0, sep='\t', engine='python')
        id_col_name = df.columns[0]
    except Exception as e:
        print(f"❌ An error occurred while reading the file: {e}")
        return

    print("Parsing identifiers...")
    temp_parsed_info = df[id_col_name].apply(parse_identifier)
    df['ParsedSequenceHeader'] = [item[0] for item in temp_parsed_info]
    df['ParsedChunkNum'] = [item[1] for item in temp_parsed_info]
    
    data_cols = df.columns.drop([id_col_name, 'ParsedSequenceHeader', 'ParsedChunkNum'])
    
    df.dropna(subset=['ParsedSequenceHeader', 'ParsedChunkNum'], inplace=True)
    df['ParsedChunkNum'] = df['ParsedChunkNum'].astype(int)
    
    grouped_by_sequence = df.groupby('ParsedSequenceHeader')
    print(f"Found {len(grouped_by_sequence)} unique sequence(s) to process.")

    master_events_list = []

    for seq_header, group_df in grouped_by_sequence:
        print(f"\r  -> Processing sequence: {seq_header[:50]:<50}", end="")
        
        for tfbs_col in data_cols:
            active_chunks = group_df[group_df[tfbs_col] == 1]['ParsedChunkNum'].sort_values()
            
            # --- MODIFIED: Ensure chunks are within the valid range based on sequence length ---
            max_chunk_num = ((TOTAL_SEQUENCE_LENGTH - WINDOW_SIZE) // STEP_SIZE) + 1
            active_chunks = active_chunks[active_chunks <= max_chunk_num]
            if active_chunks.empty:
                continue

            island_ids = (active_chunks.diff() > 1).cumsum()
            
            # --- MODIFIED BLOCK START: Gold Standard Event Splitting Logic ---
            for _, chunks_in_island in active_chunks.groupby(island_ids):
                chunk_min_num, chunk_max_num = chunks_in_island.min(), chunks_in_island.max()
                predicted_start_region = (chunk_max_num - 1) * STEP_SIZE + 1
                predicted_end_region = ((chunk_min_num - 1) * STEP_SIZE + 1) + WINDOW_SIZE - 1
                
                predicted_end_region = min(predicted_end_region, TOTAL_SEQUENCE_LENGTH)

                if predicted_start_region > predicted_end_region:
                    central_chunk_num = chunks_in_island.iloc[len(chunks_in_island)//2]
                    predicted_start_region = (central_chunk_num - 1) * STEP_SIZE + 1
                    predicted_end_region = predicted_start_region + WINDOW_SIZE - 1
                    predicted_end_region = min(predicted_end_region, TOTAL_SEQUENCE_LENGTH)

                original_pi = len(chunks_in_island) / MAX_WINDOW_OVERLAP
                remaining_pi = original_pi
                
                num_sites = int(np.ceil(original_pi)) if original_pi > 0 else 1
                total_width = float(predicted_end_region - predicted_start_region)
                site_width = total_width / num_sites if num_sites > 0 else total_width
                
                current_start = float(predicted_start_region)

                while remaining_pi > 0:
                    pi_for_this_site = min(remaining_pi, 1.0)
                    
                    site_end = current_start + site_width
                    
                    if remaining_pi <= 1.0:
                        site_end = predicted_end_region

                    if pi_for_this_site >= 1.0:
                        uncertainty = 1.0
                    else:
                        uncertainty = 250 - (240 * pi_for_this_site)

                    master_events_list.append({
                        'TFBS_Family': tfbs_col,
                        'PI': pi_for_this_site,
                        'Predicted_Start': int(current_start),
                        'Predicted_End': int(site_end),
                        'Uncertainty_Range_bp': uncertainty
                    })

                    remaining_pi -= pi_for_this_site
                    current_start = site_end
            # --- MODIFIED BLOCK END ---

    print("\n")
    if not master_events_list:
        print("❌ No TFBS events were detected in the entire file. Exiting.")
        return
        
    print(f"✅ Found a total of {len(master_events_list)} events across all sequences.")
    
    print("Calculating summary statistics...")
    all_events_df = pd.DataFrame(master_events_list)

    def high_pi_proportion(x):
        return (x >= 0.5).mean()

    def se_proportion(x):
        p = (x >= 0.5).mean()
        n = len(x)
        return np.sqrt(p * (1 - p) / n) if n > 0 else 0
    
    summary_stats = all_events_df.groupby('TFBS_Family').agg(
        Event_Count=('PI', 'size'),
        Average_PI=('PI', 'mean'),
        PI_Error=('PI', 'sem'),
        Average_Start=('Predicted_Start', 'mean'),
        Start_Error=('Predicted_Start', 'sem'),
        Average_End=('Predicted_End', 'mean'),
        End_Error=('Predicted_End', 'sem'),
        Average_Uncertainty_Range=('Uncertainty_Range_bp', 'mean'),
        Uncertainty_Range_Error=('Uncertainty_Range_bp', 'sem'),
        High_PI_Proportion=('PI', high_pi_proportion),
        High_PI_Error=('PI', se_proportion)
    ).reset_index()

    print("\n--- Summary Statistics Table ---")
    print(summary_stats.to_string())
    table_path = output_path / 'summary_statistics.csv'
    summary_stats.to_csv(table_path, index=False, float_format='%.4f')
    print("------------------------------")
    print(f"  -> ✅ Saved full summary table to: {table_path}\n")

    plt.style.use('seaborn-v0_8-whitegrid')

    # Plot 1: Total Prediction Counts
    print("📊 Generating total prediction count plot...")
    count_sorted = summary_stats.sort_values('Event_Count', ascending=False)
    plt.figure(figsize=(14, 8))
    bars = plt.bar(
        count_sorted['TFBS_Family'],
        count_sorted['Event_Count'],
        color=NATURE_PALETTE['count'],
        edgecolor='black'
    )
    plt.bar_label(bars, padding=3, fontsize=TICK_FONTSIZE-2)
    plt.title('Total Prediction Counts per TFBS Family', fontsize=TITLE_FONTSIZE, weight='bold')
    plt.xlabel('TFBS Family', fontsize=LABEL_FONTSIZE)
    plt.ylabel('Total Number of Predicted Events', fontsize=LABEL_FONTSIZE)
    plt.xticks(rotation=45, ha='right', fontsize=TICK_FONTSIZE)
    plt.yticks(fontsize=TICK_FONTSIZE)
    plt.ylim(top=count_sorted['Event_Count'].max() * 1.1)
    plt.tight_layout()
    plt.savefig(output_path / 'summary_prediction_counts_barplot.png', dpi=DPI)
    plt.close()
    count_table_path = output_path / 'prediction_counts_data.csv'
    count_sorted[['TFBS_Family', 'Event_Count']].to_csv(count_table_path, index=False)
    print(f"  -> ✅ Saved prediction count plot to: {output_path / 'summary_prediction_counts_barplot.png'}")
    print(f"  -> ✅ Saved prediction count data to: {count_table_path}")

    # Plot 2: Average Positional Independence (PI)
    print("📊 Generating Positional Independence (PI) summary plot...")
    pi_sorted = summary_stats.sort_values('Average_PI')
    plt.figure(figsize=(14, 8))
    plt.bar(
        pi_sorted['TFBS_Family'], 
        pi_sorted['Average_PI'], 
        yerr=pi_sorted['PI_Error'], 
        capsize=5, color=NATURE_PALETTE['pi'], edgecolor='black'
    )
    plt.title('Average Positional Independence (PI) of TFBS Families', fontsize=TITLE_FONTSIZE, weight='bold')
    plt.xlabel('TFBS Family', fontsize=LABEL_FONTSIZE)
    plt.ylabel('Average Positional Independence (PI)', fontsize=LABEL_FONTSIZE)
    plt.xticks(rotation=45, ha='right', fontsize=TICK_FONTSIZE)
    plt.yticks(fontsize=TICK_FONTSIZE)
    plt.ylim(bottom=0)
    plt.tight_layout()
    plt.savefig(output_path / 'summary_pi_barplot.png', dpi=DPI)
    plt.close()
    pi_table_path = output_path / 'average_pi_data.csv'
    pi_sorted[['TFBS_Family', 'Average_PI', 'PI_Error']].to_csv(pi_table_path, index=False)
    print(f"  -> ✅ Saved PI plot to: {output_path / 'summary_pi_barplot.png'}")
    print(f"  -> ✅ Saved average PI data to: {pi_table_path}")

    # Plot 3: Start & End Position Distributions
    print("📊 Generating start/end position boxplot for events with PI > 0.5 (Range < 125 bp)...")
    high_pi_events_df = all_events_df[all_events_df['PI'] > 0.5].copy()
    if high_pi_events_df.empty:
        print("  -> ⚠️ Skipping position boxplot: No events found with PI > 0.5")
    else:
        median_starts = high_pi_events_df.groupby('TFBS_Family')['Predicted_Start'].median().sort_values()
        sorted_tfbs_order = median_starts.index
        melted_df = high_pi_events_df.melt(
            id_vars=['TFBS_Family'],
            value_vars=['Predicted_Start', 'Predicted_End'],
            var_name='Position_Type',
            value_name='Position'
        )
        melted_df['Position_Type'] = melted_df['Position_Type'].str.replace('Predicted_', '').str.capitalize()
        plt.figure(figsize=(16, 9))
        sns.boxplot(
            x='TFBS_Family', y='Position', hue='Position_Type',
            data=melted_df, order=sorted_tfbs_order,
            palette=PALETTE_START_END, fliersize=2.5
        )
        plt.title('Start & End Position Distributions (PI > 0.5)', fontsize=TITLE_FONTSIZE, weight='bold')
        plt.xlabel('TFBS Family', fontsize=LABEL_FONTSIZE)
        plt.ylabel('Genomic Position (bp)', fontsize=LABEL_FONTSIZE)
        plt.xticks(rotation=45, ha='right', fontsize=TICK_FONTSIZE)
        plt.yticks(fontsize=TICK_FONTSIZE)
        plt.legend(title='Position Type', fontsize=LABEL_FONTSIZE-2)
        plt.tight_layout()
        plt.savefig(output_path / 'summary_start_end_position_boxplot.png', dpi=DPI)
        plt.close()
        
        position_summary = high_pi_events_df.groupby('TFBS_Family').agg(
            {
                'Predicted_Start': ['mean', 'std', lambda x: x.quantile(0.25), lambda x: x.quantile(0.50), lambda x: x.quantile(0.75)],
                'Predicted_End': ['mean', 'std', lambda x: x.quantile(0.25), lambda x: x.quantile(0.50), lambda x: x.quantile(0.75)]
            }
        )
        position_summary.columns = ['_'.join(col).strip() for col in position_summary.columns.values]
        position_summary = position_summary.rename(columns={
            'Predicted_Start_<lambda_0>': 'Start_q25', 'Predicted_Start_<lambda_1>': 'Start_q50_median', 'Predicted_Start_<lambda_2>': 'Start_q75',
            'Predicted_End_<lambda_0>': 'End_q25', 'Predicted_End_<lambda_1>': 'End_q50_median', 'Predicted_End_<lambda_2>': 'End_q75'
        })

        position_table_path = output_path / 'promoter_position_summary_statistics.csv'
        position_summary.to_csv(position_table_path, float_format='%.4f')
        print(f"  -> ✅ Saved start/end position boxplot to: {output_path / 'summary_start_end_position_boxplot.png'}")
        print(f"  -> ✅ Saved detailed promoter position summary statistics to: {position_table_path}")

    # Plot 4: High PI Proportion
    print("📊 Generating high-PI proportion plot...")
    high_pi_sorted = summary_stats.sort_values('High_PI_Proportion')
    plt.figure(figsize=(14, 8))
    plt.bar(
        high_pi_sorted['TFBS_Family'],
        high_pi_sorted['High_PI_Proportion'],
        yerr=high_pi_sorted['High_PI_Error'],
        capsize=5, color=NATURE_PALETTE['proportion'], edgecolor='black'
    )
    plt.title('Proportion of High-PI Events (PI > 0.5)', fontsize=TITLE_FONTSIZE, weight='bold')
    plt.xlabel('TFBS Family', fontsize=LABEL_FONTSIZE)
    plt.ylabel('Proportion of Events', fontsize=LABEL_FONTSIZE)
    plt.xticks(rotation=45, ha='right', fontsize=TICK_FONTSIZE)
    plt.yticks(fontsize=TICK_FONTSIZE)
    plt.ylim(0, max(0.3, high_pi_sorted['High_PI_Proportion'].max() * 1.1))
    plt.tight_layout()
    plt.savefig(output_path / 'summary_high_pi_barplot.png', dpi=DPI)
    plt.close()
    proportion_table_path = output_path / 'high_pi_proportion_data.csv'
    high_pi_sorted[['TFBS_Family', 'High_PI_Proportion', 'High_PI_Error']].to_csv(proportion_table_path, index=False)
    print(f"  -> ✅ Saved high-PI proportion plot to: {output_path / 'summary_high_pi_barplot.png'}")
    print(f"  -> ✅ Saved high-PI proportion data to: {proportion_table_path}")

    # Plot 5: Average Positional Uncertainty
    print("📊 Generating positional uncertainty plot...")
    uncertainty_sorted = summary_stats.sort_values('Average_Uncertainty_Range')
    plt.figure(figsize=(14, 8))
    plt.bar(
        uncertainty_sorted['TFBS_Family'],
        uncertainty_sorted['Average_Uncertainty_Range'],
        yerr=uncertainty_sorted['Uncertainty_Range_Error'],
        capsize=5, color=NATURE_PALETTE['uncertainty'], edgecolor='black'
    )
    plt.title('Average Positional Uncertainty per TFBS Family', fontsize=TITLE_FONTSIZE, weight='bold')
    plt.xlabel('TFBS Family', fontsize=LABEL_FONTSIZE)
    plt.ylabel('Average Uncertainty Range (bp)', fontsize=LABEL_FONTSIZE)
    plt.xticks(rotation=45, ha='right', fontsize=TICK_FONTSIZE)
    plt.yticks(fontsize=TICK_FONTSIZE)
    plt.ylim(bottom=0)
    plt.tight_layout()
    plt.savefig(output_path / 'summary_uncertainty_range_barplot.png', dpi=DPI)
    plt.close()
    uncertainty_table_path = output_path / 'average_uncertainty_range_data.csv'
    uncertainty_sorted[['TFBS_Family', 'Average_Uncertainty_Range', 'Uncertainty_Range_Error']].to_csv(uncertainty_table_path, index=False)
    print(f"  -> ✅ Saved positional uncertainty plot to: {output_path / 'summary_uncertainty_range_barplot.png'}")
    print(f"  -> ✅ Saved positional uncertainty data to: {uncertainty_table_path}")
    
    print("\n🎉 All analysis and plotting completed successfully.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Analyzes TFBS predictions to calculate and plot summary statistics across all binding site families.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-i', '--input_file', type=str, required=True, help="Path to the input TSV data file.")
    parser.add_argument('-o', '--output_directory', type=str, required=True, help="Path to the directory to save summary plots and table.")
    
    args = parser.parse_args()
    analyze_and_summarize_events(args.input_file, args.output_directory)
