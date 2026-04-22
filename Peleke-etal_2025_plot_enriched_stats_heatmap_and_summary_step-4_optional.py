import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import argparse
import sys
import re

# ==========================================
# CONFIGURATION
# ==========================================

#    "Generate Smart Heatmaps & Stats for DeepCIS Families."
#    "Path to enriched_stats.tsv from step 3



plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']

# --- SMART MAPPING ---
FAMILY_KEYWORD_MAPPING = {
    'AP2EREBP': ['AP2', 'EREBP', 'ERF', 'DREB', 'AP2/EREBP'],
    'bHLH':     ['bHLH', 'Helix-Loop-Helix'],
    'BBRBPC':   ['BPC','BBR', 'BBR/BPC'],
    'HB':       ['Homeobox', 'HB', 'HD-ZIP'],
    'MYB':      ['MYB', 'SANT'],
    'NAC':      ['NAC', 'NAM'],
    'WRKY':     ['WRKY'],
    'C2H2':     ['C2H2', 'Zinc finger'],
    'C2C2':     ['C2C2', 'Zinc finger'],
    'C2C2dof':  ['DOF', 'Other C4 zinc finger-type factors', 'C2C2', 'Zinc finger'],
    'G2like':   ['G2-like', 'G2like', 'GARP', 'ARR', 'GARP_ARR']
}

def get_keywords_from_family(family_name):
    """Returns search keywords. Checks manual map first, then falls back to splitting."""
    s_name = str(family_name).strip()
    
    if s_name in FAMILY_KEYWORD_MAPPING:
        return FAMILY_KEYWORD_MAPPING[s_name]
    
    tokens = re.split(r'[_\W]+', s_name)
    return [t for t in tokens if t and len(t) > 1]

def filter_jaspar_by_target(df, target_family):
    keywords = get_keywords_from_family(target_family)
    if not keywords:
        return pd.DataFrame()
    
    pattern = '|'.join(map(re.escape, keywords))
    
    mask_family = df['JASPAR_tf_family'].fillna('').astype(str).str.contains(pattern, case=False, regex=True)
    mask_class  = df['tf_class'].fillna('').astype(str).str.contains(pattern, case=False, regex=True)
    
    return df[mask_family | mask_class].copy()

def plot_heatmap(df_subset, target_family, output_path, global_type_color_map):
    # 1. PREPARE LABELS & SORTING
    df_subset['safe_tf_name'] = df_subset['tf_name'].fillna('Unknown')
    df_subset['Label'] = df_subset['TF_id'] + " (" + df_subset['safe_tf_name'].astype(str) + ")"
    
    df_subset.sort_values(by=['data_type', 'safe_tf_name', 'TF_id'], ascending=[True, True, True], inplace=True)
    
    ordered_labels = df_subset['Label'].drop_duplicates().tolist()

    # 2. PIVOT & REINDEX
    heatmap_data = df_subset.pivot_table(index='Label', columns='Predicted_Family', values='Recovery_Rate')
    heatmap_data = heatmap_data.fillna(0)
    heatmap_data = heatmap_data.reindex(ordered_labels)
    heatmap_data = heatmap_data.sort_index(axis=1)

    if heatmap_data.empty:
        return

    # 3. ROW COLORS
    label_to_datatype = df_subset.drop_duplicates(subset='Label').set_index('Label')['data_type']
    row_datatypes = label_to_datatype.loc[heatmap_data.index]
    row_colors = row_datatypes.map(lambda x: global_type_color_map.get(x, '#d3d3d3'))

    # 4. DYNAMIC HEIGHT
    num_rows = len(heatmap_data)
    calculated_height = max(6, num_rows * 0.25)

    # 5. PLOT CONFIGURATION
    plt.figure(figsize=(16, calculated_height))
    cbar_left_pos = 1.02
    
    g = sns.clustermap(
        heatmap_data,
        row_colors=row_colors,
        row_cluster=False,
        col_cluster=False,
        cmap="viridis",
        figsize=(16, calculated_height),
        yticklabels=True,
        xticklabels=True,
        dendrogram_ratio=(0.01, 0.01),
        cbar_pos=(cbar_left_pos, 0.05, 0.02, 0.20) 
    )

    # 6. FORMATTING
    y_fontsize = 10 if num_rows < 50 else 8
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=y_fontsize)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize=10, rotation=45, ha='right')
    
    g.ax_heatmap.set_xlabel("Predicted Family (DeepCIS)", fontsize=12, fontweight='bold')
    g.ax_heatmap.set_ylabel("JASPAR TF ID (Name)", fontsize=12, fontweight='bold')
    g.ax_heatmap.set_title(f"Recall Rates for {target_family} Family", fontsize=14, pad=20)

    # 7. LEGEND
    legend_patches = [mpatches.Patch(color=color, label=str(label)) 
                      for label, color in global_type_color_map.items()]
    
    g.ax_heatmap.legend(
        handles=legend_patches, 
        title="Data Type", 
        loc="lower left",            
        bbox_to_anchor=(1.15, 0.0), 
        borderaxespad=0.,
        frameon=False
    )

    if hasattr(g, 'cbar_ax'):
        g.cbar_ax.set_ylabel('Recall (Recovery Rate)', rotation=270, labelpad=20, fontsize=11)

    # 8. SAVE
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close() 

def calculate_stats(df_target, family_name):
    """
    Calculates Avg, Min, Max for Sensitivity (Recovery_Rate), Precision, and F1.
    """
    stats = {
        'Family': family_name,
        'Num_Motifs': len(df_target)
    }

    # Map output names to expected column names in input file
    # Ensure these column names match your input file!
    metrics_map = {
        'Sensitivity': 'Recovery_Rate',
        'Precision': 'Precision',  # Will check if this exists
        'F1': 'F1'                 # Will check if this exists
    }

    for metric_label, col_name in metrics_map.items():
        if col_name in df_target.columns:
            # Drop NaNs just for calculation to be safe
            vals = df_target[col_name].dropna()
            if not vals.empty:
                stats[f'Avg_{metric_label}'] = round(vals.mean(), 4)
                stats[f'Min_{metric_label}'] = round(vals.min(), 4)
                stats[f'Max_{metric_label}'] = round(vals.max(), 4)
            else:
                stats[f'Avg_{metric_label}'] = 0.0
                stats[f'Min_{metric_label}'] = 0.0
                stats[f'Max_{metric_label}'] = 0.0
        else:
            # If column is missing, set to None or 0
            stats[f'Avg_{metric_label}'] = "N/A"
            stats[f'Min_{metric_label}'] = "N/A"
            stats[f'Max_{metric_label}'] = "N/A"

    return stats

def main():
    parser = argparse.ArgumentParser(description="Generate Smart Heatmaps & Stats for DeepCIS Families.")
    parser.add_argument("-i", "--input", required=True, help="Path to enriched_stats.tsv")
    parser.add_argument("-o", "--outdir", required=True, help="Directory to save plots and stats")
    
    args = parser.parse_args()

    if not os.path.exists(args.input):
        sys.exit(f"Error: Input file {args.input} not found.")
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    print(f"[INFO] Loading data from {args.input}...")
    df = pd.read_csv(args.input, sep='\t')
    df['data_type'] = df['data_type'].fillna('Unknown')

    # GLOBAL COLORS
    all_data_types = sorted(df['data_type'].unique())
    palette = sns.color_palette("husl", len(all_data_types))
    global_type_color_map = dict(zip(all_data_types, palette))

    target_families = df['Predicted_Family'].dropna().unique()
    print(f"[INFO] Processing {len(target_families)} target families...")

    # List to collect stats for all families
    all_stats = []

    for target in target_families:
        safe_name = re.sub(r'[^\w\-_]', '_', str(target))
        out_file = os.path.join(args.outdir, f"heatmap_{safe_name}.png")
        
        df_target = filter_jaspar_by_target(df, target)
        if df_target.empty:
            continue
        
        # --- 1. CALCULATE STATS ---
        family_stats = calculate_stats(df_target, target)
        all_stats.append(family_stats)

        # --- 2. PLOT HEATMAP ---
        num_rows = len(df_target)
        print(f" > Processing {target}: {num_rows} motifs found.")
        try:
            plot_heatmap(df_target, target, out_file, global_type_color_map)
        except Exception as e:
            print(f"   [ERROR] Failed to plot {target}: {e}")

    # --- 3. SAVE AND PRINT SUMMARY STATS ---
    if all_stats:
        stats_df = pd.DataFrame(all_stats)
        
        # Define column order for cleaner output
        cols_order = ['Family', 'Num_Motifs', 
                      'Avg_Sensitivity', 'Min_Sensitivity', 'Max_Sensitivity',
                      'Avg_Precision', 'Min_Precision', 'Max_Precision',
                      'Avg_F1', 'Min_F1', 'Max_F1']
        
        # Filter cols_order to only include columns that actually exist in stats_df
        final_cols = [c for c in cols_order if c in stats_df.columns]
        stats_df = stats_df[final_cols]

        # Save to CSV
        stats_path = os.path.join(args.outdir, "family_summary_stats.csv")
        stats_df.to_csv(stats_path, index=False)
        
        print("\n" + "="*60)
        print(" SUMMARY STATISTICS")
        print("="*60)
        # Convert to string to prevent truncation in console view
        print(stats_df.to_string(index=False))
        print("="*60)
        print(f"[SUCCESS] Statistics saved to: {stats_path}")
    
    print("[SUCCESS] All tasks completed.")

if __name__ == "__main__":
    main()
