#!/usr/bin/env python3
import argparse
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from adjustText import adjust_text

# --- Constants ---
TOTAL_SEQUENCE_LENGTH = 1500
WINDOW_SIZE = 250
STEP_SIZE = 10
MAX_WINDOW_OVERLAP = WINDOW_SIZE // STEP_SIZE

def parse_identifier(identifier_str):
    """Parses an identifier string to extract the base sequence identifier and chunk number."""
    match = re.search(r'(.*)_chunk(\d+)', str(identifier_str))
    if match: return match.group(1), int(match.group(2))
    return None, None

def sanitize_filename(name):
    """Removes characters from a string that are invalid for filenames."""
    return re.sub(r'[^a-zA-Z0-9_.-]', '_', name)

def analyze_and_pinpoint_events(input_file, output_dir):
    """
    Main function to identify and characterize putative TFBS event regions.
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

    print("Generating extended color map...")
    all_tfbs_families = data_cols.tolist()
    num_families = len(all_tfbs_families)
    colors_tab20 = plt.cm.tab20(np.linspace(0, 1, 20))
    colors_tab20b = plt.cm.tab20b(np.linspace(0, 1, 20))
    colors_set3 = plt.cm.Set3(np.linspace(0, 1, 12))
    extended_palette = np.vstack((colors_tab20, colors_tab20b, colors_set3))
    fixed_colors = [extended_palette[i % len(extended_palette)] for i in range(num_families)]
    global_color_map = dict(zip(all_tfbs_families, fixed_colors))

    df.dropna(subset=['ParsedSequenceHeader', 'ParsedChunkNum'], inplace=True)
    df['ParsedChunkNum'] = df['ParsedChunkNum'].astype(int)

    grouped_by_sequence = df.groupby('ParsedSequenceHeader')
    print(f"Found {len(grouped_by_sequence)} unique sequence(s) to process.")

    for seq_header, group_df in grouped_by_sequence:
        print(f"\nProcessing sequence: {seq_header}...")
        all_events = []

        for tfbs_col in data_cols:
            active_chunks = group_df[group_df[tfbs_col] == 1]['ParsedChunkNum'].sort_values()
            
            max_chunk_num = ((TOTAL_SEQUENCE_LENGTH - WINDOW_SIZE) // STEP_SIZE) + 1
            active_chunks = active_chunks[active_chunks <= max_chunk_num]
            if active_chunks.empty:
                continue

            island_ids = (active_chunks.diff() > 1).cumsum()

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

                    # --- MODIFIED BLOCK START: Corrected Position & Uncertainty Logic ---
                    if pi_for_this_site >= 1.0:
                        uncertainty = 1.0
                        # High PI: Pinpoint the position near the start of the site.
                        predicted_position_for_site = current_start + 10
                    else:
                        uncertainty = 250 - (240 * pi_for_this_site)
                        # Low PI: Anchor the event to the start of the site.
                        predicted_position_for_site = current_start
                    # --- MODIFIED BLOCK END ---

                    all_events.append({
                        'TFBS_Family': tfbs_col,
                        'Predicted_Position': predicted_position_for_site,
                        'PI': pi_for_this_site,
                        'Predicted_Start': int(current_start),
                        'Predicted_End': int(site_end),
                        'Uncertainty_Range_bp': uncertainty
                    })

                    remaining_pi -= pi_for_this_site
                    current_start = site_end

        if not all_events:
            print("  -> No TFBS events detected for this sequence.")
            continue

        events_df = pd.DataFrame(all_events)
        events_df.sort_values(by=['Predicted_Start'], inplace=True)
        
        events_df['TFBS_Number'] = events_df.groupby('TFBS_Family').cumcount() + 1
        
        safe_seq_header = sanitize_filename(seq_header)
        table_filename = output_path / f"{safe_seq_header}_predicted_events.csv"

        cols_to_save = ['TFBS_Family', 'TFBS_Number', 'Predicted_Position', 'PI',
                        'Predicted_Start', 'Predicted_End', 'Uncertainty_Range_bp']
        events_df[cols_to_save].to_csv(table_filename, index=False, float_format='%.4f')
        print(f"  -> ✅ Found {len(events_df)} events. Saved summary table to: {table_filename}")

        plt.style.use('seaborn-v0_8-whitegrid')
        fig, ax = plt.subplots(figsize=(12, 2.5))

        ax.grid(which='major', axis='y', linestyle=':', linewidth='0.5', color='lightgray')

        plot_df = events_df[events_df['PI'] > 0].copy()
        
        texts = []
        for _, event in plot_df.iterrows():
            # The triangle is now centered on the corrected, start-anchored position
            center_pos = event['Predicted_Position']
            triangle_width = event['Uncertainty_Range_bp'] + 10
            
            x_start = center_pos - (triangle_width / 2)
            x_end = center_pos + (triangle_width / 2)
            
            x_peak = center_pos
            y_peak = event['PI']

            x_coords = [x_start, x_peak, x_end]
            y_coords = [0, y_peak, 0]
            
            ax.fill_between(
                x_coords,
                y_coords,
                color=global_color_map[event['TFBS_Family']],
                alpha=0.7,
                zorder=10
            )
            
            label_text = f"{event['TFBS_Family']}_{event['TFBS_Number']}"
            texts.append(ax.text(
                x=x_peak,
                y=y_peak + 0.02, # Start text just above the peak
                s=label_text,
                ha='center',
                fontsize=6
            ))
        
        adjust_text(
            texts,
            ax=ax,
            arrowprops=dict(
                arrowstyle='-', 
                color='gray',
                lw=0.5
            )
        )

        tss_position = 1000
        ax.axvline(x=tss_position, color='crimson', linestyle='--', linewidth=1.2, zorder=5)
        
        ax.set_title(f"Predicted TFBS Events for {seq_header}", fontsize=9, weight='bold', pad=40)
        ax.set_xlabel("Genomic Position (bp)", fontsize=8)

        ax.set_ylabel("Positional Independence (PI)", fontsize=8)
        ax.set_ylim(0, 1.05)
        ax.tick_params(axis='y', which='major', labelsize=7)
        ax.text(tss_position, ax.get_ylim()[1] + 0.1, 'TSS', # Place TSS label above y-axis
                color='crimson', ha='center', va='bottom', fontsize=8, weight='bold')

        ax2 = ax.twinx()
        ax2.set_ylim(250 + 10, 1 + 10)
        ax2.set_ylabel("Binding Precision (bp)", fontsize=8, color='royalblue', fontweight='bold')
        ax2.tick_params(axis='y', colors='royalblue', labelsize=7)

        major_ticks = np.arange(0, TOTAL_SEQUENCE_LENGTH + 1, 50)
        ax.set_xticks(major_ticks)
        ax.set_xticks([], minor=True)
        
        ax.tick_params(axis='x', which='major', labelsize=6, length=4, direction='out', top=False, rotation=90)

        ax.set_xlim(0, TOTAL_SEQUENCE_LENGTH + 1)

        active_tfbs_families = plot_df['TFBS_Family'].unique()
        handles = [plt.Rectangle((0,0),1,1, color=global_color_map[fam], alpha=0.7) for fam in active_tfbs_families]
        
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
        
        ax.legend(handles, active_tfbs_families,
                    loc='center left',
                    bbox_to_anchor=(1.05, 0.5),
                    fontsize=6,
                    title='TFBS Families',
                    title_fontsize=7,
                    ncol=2
                   )
        
        plt.tight_layout(rect=[0, 0, 1, 0.95])

        plot_filename = output_path / f"{safe_seq_header}_events_plot_final.png"
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"  -> ✅ Saved final plot to: {plot_filename}")

    print("\n🎉 All sequences processed successfully.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Detects and characterizes putative TFBS binding event regions from overlapping window predictions.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-i', '--input_file', type=str, required=True, help="Path to the input data file.")
    parser.add_argument('-o', '--output_directory', type=str, required=True, help="Path to the directory for output files.")
    
    try:
        import adjustText
    except ImportError:
        print("\n---")
        print("Warning: The 'adjustText' library is not installed.")
        print("This script requires it for non-overlapping label placement.")
        print("Please install it by running: pip install adjustText")
        print("---\n")
        exit()
        
    args = parser.parse_args()
    analyze_and_pinpoint_events(args.input_file, args.output_directory)
