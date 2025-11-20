import pandas as pd
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import io
from adjustText import adjust_text # Import for non-overlapping labels

# --- 1. Configuration ---

# --- File Paths based on execution location ---
# The script is executed from: ~/Desktop/deepCIS_calc/workflows
DATA_FILE_PATH = '../studies/deepCIS/deepCistome/EVALUATION/TFBS_performance_metrics_complete_sz25.csv'
OUTPUT_DIR = 'new_linear_regression_plots' # Changed output directory
OUTPUT_FILENAME = os.path.join(OUTPUT_DIR, 'combined_regression_plot.png')

# Define the pairs of *DISPLAY LABELS* for regression
REGRESSION_PAIRS = [
    ('TF-WOI', 'sensitivity'),
    ('TF-WOI', 'MCC'),
    ('IPM predictability', 'sensitivity'),
    ('IPM predictability', 'MCC'),
    ('TF-CoO', 'sensitivity'),
    ('TF-CoO', 'MCC')
]

# *** NEW: Map display labels to the actual column names in the CSV ***
COLUMN_NAME_MAP = {
    # Display Label : Technical Column Name
    'TF-WOI': 'aver_TF_PI',
    'IPM predictability': 'aver_TF_IPM_pred',
    'TF-CoO': 'w_aver_TF_CoOcc',
    'sensitivity': 'sensitivity', # This one is the same
    'MCC': 'MCC'               # This one is the same
}

# --- 2. Setup Plotting Style ---

# Use very large fonts for a compact, readable plot
sns.set_style("ticks")
plt.rcParams.update({
    'font.size': 20,
    'axes.labelsize': 16, # Adjusted axis label size
    'xtick.labelsize': 14, # Adjusted tick label size
    'ytick.labelsize': 14, # Adjusted tick label size
    'figure.titlesize': 26,
    'axes.titleweight': 'bold',
})

# --- 3. Main Analysis Function ---

def perform_regression_and_plot():
    """
    Loads data, performs linear regression, and
    generates a single combined plot for all specified contrasts.
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

    print(f"\n--- Starting Linear Regression ---") # Updated print message

    # --- Create a single figure with a 3x2 grid of subplots ---
    # Increased height to give plots more space
    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(10, 15))
    axes = axes.flatten() # Flatten the 2D array of axes for easy iteration

    # --- Loop Through Each Regression Task and plot on a subplot ---
    # We now loop through the *labels*
    for i, (y_label, x_label) in enumerate(REGRESSION_PAIRS):
        ax = axes[i] # Get the current subplot
        
        # *** NEW: Get the technical column names from the map ***
        try:
            y_col = COLUMN_NAME_MAP[y_label]
            x_col = COLUMN_NAME_MAP[x_label]
        except KeyError as e:
            print(f"❌ Error: Label '{e}' not found in COLUMN_NAME_MAP. Skipping plot.")
            ax.text(0.5, 0.5, f"Error: No data for '{e}'", ha='center', va='center', fontsize=14, color='red')
            continue

        print(f"\n📊 Plotting on subplot {i+1}: '{y_label}' (Y) vs. '{x_label}' (X)")

        # *** MODIFIED: Use technical names (y_col, x_col) to find data ***
        analysis_df = df.dropna(subset=[y_col, x_col, 'TF']).copy()

        if analysis_df.empty:
            ax.text(0.5, 0.5, 'No valid data', ha='center', va='center', fontsize=18)
            # *** MODIFIED: Use display labels (y_label, x_label) for the empty plot title ***
            ax.set_title(f"{y_label}\nvs.\n{x_label}", fontsize=18)
            sns.despine()
            continue

        # --- Perform Linear Regression ---
        # *** MODIFIED: Use technical names (y_col, x_col) for the formula ***
        formula = f"{y_col} ~ {x_col}"
        
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

        # --- Plotting on the current subplot ---
        # *** MODIFIED: Use technical names (y_col, x_col) to plot data points ***
        ax.scatter(inliers_df[x_col], inliers_df[y_col], color='darkgrey', s=60, alpha=0.3)
        ax.scatter(outliers_df[x_col], outliers_df[y_col], color='orange', s=70)

        # Generate the regression line
        # *** MODIFIED: Use technical names (y_col, x_col) for prediction ***
        x_curve = pd.Series(sorted(analysis_df[x_col].unique()))
        y_curve = model.predict(pd.DataFrame({x_col: x_curve}))
        ax.plot(x_curve, y_curve, color='black', linewidth=3) # Color was 'black' in your script

        # Collect text labels for automatic adjustment
        texts = []
        for _, row in outliers_df.iterrows():
            # *** MODIFIED: Use technical names (y_col, x_col) to position text ***
            texts.append(ax.text(row[x_col], row[y_col], f" {row['TF']}", fontsize=12))

        # Automatically adjust text labels to prevent overlap
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

        # --- Final subplot styling ---
        # Move stats into the title
        stats_text = f'$R^2={r_squared:.2f}, p={p_value:.1e}$'
        ax.set_title(f'{stats_text}', fontsize=12)

        sns.despine()
        # *** MODIFIED: Use display labels (y_label, x_label) for axes ***
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_xlim(0, 1)

    # --- Final Figure Adjustments ---
    plt.tight_layout(pad=3.0, h_pad=4.0) # Add padding between plots
    plt.savefig(OUTPUT_FILENAME, dpi=300, bbox_inches='tight')
    plt.close(fig)

    print(f"\n✅ Analysis complete. Combined plot saved to '{OUTPUT_FILENAME}'")

# --- 4. Execute the Script ---
if __name__ == '__main__':
    perform_regression_and_plot()
