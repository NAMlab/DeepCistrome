import pandas as pd
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import io
from adjustText import adjust_text # Import for non-overlapping labels

# --- 1. Configuration ---

# Set the degree for the polynomial regression
POLYNOMIAL_DEGREE = 2
GLOBAL_FONT_SIZE = 14 # Set a single font size for all text

# --- File Paths based on execution location ---
# The script is executed from: ~/Desktop/deepCIS_calc/workflows
DATA_FILE_PATH = '../studies/deepCIS/deepCistome/EVALUATION/TFBS_performance_metrics_complete_sz25.csv'
OUTPUT_DIR = f'new_regression_plots_deg{POLYNOMIAL_DEGREE}'
OUTPUT_FILENAME = os.path.join(OUTPUT_DIR, 'combined_regression_plot.png')
# Filename for the MCC-only combined plot
MCC_OUTPUT_FILENAME = os.path.join(OUTPUT_DIR, 'combined_MCC_plot.png')


# Define the pairs of *DISPLAY LABELS* for regression
REGRESSION_PAIRS = [
    ('TF-WOI', 'sensitivity'),
    ('TF-WOI', 'MCC'),
    ('IPM predictability', 'sensitivity'),
    ('IPM predictability', 'MCC'),
    ('TF-CoO', 'sensitivity'),
    ('TF-CoO', 'MCC')
]

# Map display labels to the actual column names in the CSV
COLUMN_NAME_MAP = {
    # Display Label : Technical Column Name
    'TF-WOI': 'aver_TF_PI',
    'IPM predictability': 'aver_TF_IPM_pred',
    'TF-CoO': 'w_aver_TF_CoOcc',
    'sensitivity': 'sensitivity', # This one is the same
    'MCC': 'MCC'               # This one is the same
}

# --- 2. Setup Plotting Style ---

# Use a single font size for all plot elements
sns.set_style("ticks")
plt.rcParams.update({
    'font.size': GLOBAL_FONT_SIZE,
    'axes.labelsize': GLOBAL_FONT_SIZE,
    'xtick.labelsize': GLOBAL_FONT_SIZE,
    'ytick.labelsize': GLOBAL_FONT_SIZE,
    'figure.titlesize': GLOBAL_FONT_SIZE,
    'axes.titleweight': 'bold',
})

# --- 3. Plotting Helper Function ---

def create_plot_on_axis(ax, analysis_df, inliers_df, outliers_df, model, 
                        x_col, y_col, x_label, y_label, r_squared, p_value):
    """
    Helper function to draw the complete scatter/regression plot on a given axis.
    """
    # --- Plotting Data Points ---
    ax.scatter(inliers_df[x_col], inliers_df[y_col], color='darkgrey', s=60, alpha=0.3)
    ax.scatter(outliers_df[x_col], outliers_df[y_col], color='orange', s=70)

    # --- Plotting Regression Line ---
    x_curve = pd.Series(sorted(analysis_df[x_col].unique()))
    y_curve = model.predict(pd.DataFrame({x_col: x_curve}))
    ax.plot(x_curve, y_curve, color='black', linewidth=2)

    # --- Adjusting Text Labels ---
    texts = []
    for _, row in outliers_df.iterrows():
        texts.append(ax.text(row[x_col], row[y_col], f" {row['TF']}", fontsize=GLOBAL_FONT_SIZE))
    
    if texts: # Only run adjust_text if there are outliers to label
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

    # --- Final Subplot Styling ---
    # Stats text is now bold using \mathbf
    stats_text = rf'$\mathbf{{R^2={r_squared:.2f}, p={p_value:.1e}}}$'
    ax.set_title(f'{stats_text}', fontsize=GLOBAL_FONT_SIZE)

    sns.despine()
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    
    # *** THIS IS THE KEY CHANGE ***
    # Set x-axis limits conditionally
    if x_label == 'MCC':
        # For MCC, set limit from 0 to the max value + 5% buffer
        max_val = analysis_df[x_col].max()
        ax.set_xlim(0, max_val * 1.05)
    else:
        # For other plots (sensitivity), keep the 0-1 limit
        ax.set_xlim(0, 1)


# --- 4. Main Analysis Function ---

def perform_regression_and_plot():
    """
    Loads data, performs polynomial regression, and generates:
    1. A single combined plot for all contrasts.
    2. A combined plot for only MCC contrasts.
    3. Individual plots for each MCC contrast.
    """
    # --- Load and Prepare Data ---
    try:
        df = pd.read_csv(DATA_FILE_PATH)
    except FileNotFoundError:
        print(f"❌ Error: Data file not found at '{DATA_FILE_PATH}'")
        print("Please ensure the script is run from '~/Desktop/deepCIS_calc/workflows'")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error: Could not parse the data file. Error: {e}")
        sys.exit(1)

    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        print(f"📂 Created directory for plots: '{OUTPUT_DIR}'")

    print(f"\n--- Starting Degree-{POLYNOMIAL_DEGREE} Polynomial Regression (Raw Polynomials) ---")

    # --- Create Figure 1: Full 3x2 combined plot ---
    fig_all, axes_all = plt.subplots(nrows=3, ncols=2, figsize=(10, 15))
    axes_all = axes_all.flatten() # Flatten for easy iteration

    # --- Create Figure 2: MCC-only 3x1 combined plot ---
    fig_mcc, axes_mcc = plt.subplots(nrows=3, ncols=1, figsize=(5, 15))
    mcc_plot_index = 0 # Counter for the MCC subplots

    # --- Loop Through Each Regression Task ---
    for i, (y_label, x_label) in enumerate(REGRESSION_PAIRS):
        
        # --- Get technical column names ---
        try:
            y_col = COLUMN_NAME_MAP[y_label]
            x_col = COLUMN_NAME_MAP[x_label]
        except KeyError as e:
            print(f"❌ Error: Label '{e}' not found in COLUMN_NAME_MAP. Skipping plot.")
            ax = axes_all[i]
            ax.text(0.5, 0.5, f"Error: No data for '{e}'", ha='center', va='center', fontsize=GLOBAL_FONT_SIZE, color='black')
            continue
            
        print(f"\n📊 Plotting: '{y_label}' (Y) vs. '{x_label}' (X)")

        # --- Prepare data for this specific plot ---
        analysis_df = df.dropna(subset=[y_col, x_col, 'TF']).copy()

        # Handle empty data case
        if analysis_df.empty:
            print("   ...No valid data for this pair. Skipping.")
            ax = axes_all[i] # Get the current subplot for the main figure
            ax.text(0.5, 0.5, 'No valid data', ha='center', va='center', fontsize=GLOBAL_FONT_SIZE)
            ax.set_title(f"{y_label}\nvs.\n{x_label}", fontsize=GLOBAL_FONT_SIZE)
            sns.despine()
            
            # Also add empty text to MCC plot if relevant
            if x_label == 'MCC':
                # Check if axes_mcc is a single object or an array
                if isinstance(axes_mcc, plt.Axes):
                    ax_mcc = axes_mcc # Should not happen with nrows=3
                else:
                    ax_mcc = axes_mcc[mcc_plot_index]
                    
                ax_mcc.text(0.5, 0.5, 'No valid data', ha='center', va='center', fontsize=GLOBAL_FONT_SIZE)
                ax_mcc.set_title(f"{y_label}\nvs.\n{x_label}", fontsize=GLOBAL_FONT_SIZE)
                sns.despine()
                mcc_plot_index += 1
            continue

        # --- Perform Regression ---
        poly_terms = ' + '.join([f"I({x_col}**{deg})" for deg in range(2, POLYNOMIAL_DEGREE + 1)])
        formula = f"{y_col} ~ {x_col} + {poly_terms}"
        
        model = smf.ols(formula, data=analysis_df).fit()
        r_squared = model.rsquared
        p_value = model.f_pvalue

        # --- Outlier Detection ---
        residuals = model.resid
        q1, q3 = residuals.quantile(0.25), residuals.quantile(0.75)
        iqr = q3 - q1
        lower_bound = q1 - 1 * iqr
        upper_bound = q3 + 1 * iqr
        
        is_outlier = (residuals < lower_bound) | (residuals > upper_bound)
        outliers_df = analysis_df[is_outlier]
        inliers_df = analysis_df[~is_outlier]

        # --- Plot 1: Add to the main 3x2 combined figure ---
        ax_main = axes_all[i]
        create_plot_on_axis(ax_main, analysis_df, inliers_df, outliers_df, model,
                            x_col, y_col, x_label, y_label, r_squared, p_value)

        # --- Check if this is an MCC plot to create extra plots ---
        if x_label == 'MCC':
            
            # --- Plot 2: Create a new, individual plot figure ---
            print(f"   ...Creating individual plot for MCC.")
            fig_indiv, ax_indiv = plt.subplots(figsize=(6, 5)) # Create a new figure
            create_plot_on_axis(ax_indiv, analysis_df, inliers_df, outliers_df, model,
                                x_col, y_col, x_label, y_label, r_squared, p_value)
            
            # Save the individual figure
            indiv_filename = os.path.join(OUTPUT_DIR, f'individual_{y_label}_vs_{x_label}.png')
            fig_indiv.savefig(indiv_filename, dpi=300, bbox_inches='tight')
            plt.close(fig_indiv) # Close it so it doesn't stay in memory

            # --- Plot 3: Add to the MCC-only 3x1 combined figure ---
            # Check if axes_mcc is a single object or an array
            if isinstance(axes_mcc, plt.Axes):
                ax_mcc = axes_mcc # Should not happen with nrows=3
            else:
                ax_mcc = axes_mcc[mcc_plot_index]
                
            create_plot_on_axis(ax_mcc, analysis_df, inliers_df, outliers_df, model,
                                x_col, y_col, x_label, y_label, r_squared, p_value)
            mcc_plot_index += 1 # Move to the next subplot in the MCC figure

    # --- Final Figure Adjustments and Saving ---
    
    # Save the main 3x2 combined plot
    print(f"\n✅ Saving combined plot to '{OUTPUT_FILENAME}'")
    fig_all.tight_layout(pad=3.0, h_pad=4.0)
    fig_all.savefig(OUTPUT_FILENAME, dpi=300, bbox_inches='tight')
    plt.close(fig_all)

    # Save the MCC-only 3x1 combined plot
    # Add a check in case there were no MCC plots at all
    if mcc_plot_index > 0:
        print(f"✅ Saving MCC-only combined plot to '{MCC_OUTPUT_FILENAME}'")
        fig_mcc.tight_layout(pad=3.0, h_pad=4.0)
        fig_mcc.savefig(MCC_OUTPUT_FILENAME, dpi=300, bbox_inches='tight')
    plt.close(fig_mcc)

    print("\n🎉 Analysis complete.")

# --- 5. Execute the Script ---
if __name__ == '__main__':
    perform_regression_and_plot()
