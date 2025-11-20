#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import os
from scipy.stats import ttest_ind, levene
from statsmodels.stats.multitest import fdrcorrection
import statsmodels.api as sm
import statsmodels.formula.api as smf

def bootstrap_median_diff(ctrl_reps, var_reps, n_bootstraps=10000):
    """
    Performs a bootstrap analysis on the difference of medians.
    Returns the lower and upper bounds of the 95% confidence interval.
    """
    boot_diffs = np.zeros(n_bootstraps)
    for i in range(n_bootstraps):
        boot_ctrl = np.random.choice(ctrl_reps, size=len(ctrl_reps), replace=True)
        boot_var = np.random.choice(var_reps, size=len(var_reps), replace=True)
        boot_diffs[i] = np.median(boot_var) - np.median(boot_ctrl)
    
    lower_bound = np.percentile(boot_diffs, 2.5)
    upper_bound = np.percentile(boot_diffs, 97.5)
    
    return lower_bound, upper_bound

def analyze_data_comprehensive(input_file, output_file, separator):
    """
    Transforms data and runs classical, bootstrap, and regression statistical analyses.
    """
    # --- 1. Load and Prepare Data ---
    print(f"Reading data from: {input_file}")
    df = pd.read_csv(input_file, sep=separator)
    
    df['base_id'] = df['id'].str.replace(r'(-ctrl\.0|-var\.0)$', '', regex=True)
    df['condition_type'] = np.where(df['id'].str.contains('-ctrl'), 'ctrl', 'var')

    # --- 2. Reshape Data and Calculate Descriptive Statistics ---
    print("Reshaping data and calculating descriptive statistics...")
    wide_df = df.pivot_table(
        index='base_id',
        columns=['condition_type', 'rep'],
        values='enrichment'
    )
    
    wide_df.columns = ['_'.join(map(str, col)) for col in wide_df.columns.values]
    wide_df = wide_df.reset_index()

    wide_df = wide_df.dropna(subset=['ctrl_1', 'ctrl_2', 'var_1', 'var_2'])

    wide_df['ctrl_delta'] = wide_df['ctrl_1'] - wide_df['ctrl_2']
    wide_df['ctrl_var'] = wide_df[['ctrl_1', 'ctrl_2']].var(axis=1, ddof=1)
    wide_df['var_delta'] = wide_df['var_1'] - wide_df['var_2']
    wide_df['var_var'] = wide_df[['var_1', 'var_2']].var(axis=1, ddof=1)
    wide_df['ctrl_median'] = wide_df[['ctrl_1', 'ctrl_2']].median(axis=1)
    wide_df['var_median'] = wide_df[['var_1', 'var_2']].median(axis=1)
    wide_df['delta_med'] = wide_df['var_median'] - wide_df['ctrl_median']
    wide_df['var_med'] = wide_df[['ctrl_median', 'var_median']].var(axis=1, ddof=1)

    # --- 3. Perform Classical Statistical Tests ---
    print("Performing Levene and T-tests...")
    levene_pvals = []
    ttest_pvals = []
    
    for index, row in wide_df.iterrows():
        ctrl_reps = np.array([row['ctrl_1'], row['ctrl_2']])
        var_reps = np.array([row['var_1'], row['var_2']])
        
        _, p_levene = levene(ctrl_reps, var_reps)
        levene_pvals.append(p_levene)
        
        _, p_ttest = ttest_ind(ctrl_reps, var_reps, equal_var=True, nan_policy='omit')
        ttest_pvals.append(p_ttest)
        
    wide_df['p_levene'] = levene_pvals
    wide_df['p_ttest'] = ttest_pvals
    
    # FDR Correction
    if 'p_ttest' in wide_df.columns and not wide_df['p_ttest'].empty:
        _, fdr_values = fdrcorrection(wide_df['p_ttest'].dropna(), alpha=0.05)
        wide_df.loc[wide_df['p_ttest'].notna(), 'fdr'] = fdr_values

    # --- 4. Perform Bootstrap Analysis ---
    print("Performing bootstrap analysis for confidence intervals...")
    boot_results = []
    for index, row in wide_df.iterrows():
        ctrl_reps = np.array([row['ctrl_1'], row['ctrl_2']])
        var_reps = np.array([row['var_1'], row['var_2']])
        
        lower, upper = bootstrap_median_diff(ctrl_reps, var_reps)
        boot_results.append({'lower_ci': lower, 'upper_ci': upper})
        
    boot_df = pd.DataFrame(boot_results, index=wide_df.index)
    wide_df = pd.concat([wide_df, boot_df], axis=1)

    wide_df['is_significant_bootstrap'] = ~((wide_df['lower_ci'] <= 0) & (wide_df['upper_ci'] >= 0))

    # --- 5. Perform Regression Analysis ---
    print("Performing regression analysis and determining confidence intervals...")

    # Regression 1: ctrl_1 vs var_1
    model1 = smf.ols('var_1 ~ ctrl_1', data=wide_df).fit()
    pred1 = model1.get_prediction(wide_df).summary_frame(alpha=0.05)
    wide_df['ci_var1_vs_ctrl1_lower'] = pred1['mean_ci_lower']
    wide_df['ci_var1_vs_ctrl1_upper'] = pred1['mean_ci_upper']
    wide_df['within_ci_var1'] = (wide_df['var_1'] >= wide_df['ci_var1_vs_ctrl1_lower']) & \
                                (wide_df['var_1'] <= wide_df['ci_var1_vs_ctrl1_upper'])

    # Regression 2: ctrl_2 vs var_2
    model2 = smf.ols('var_2 ~ ctrl_2', data=wide_df).fit()
    pred2 = model2.get_prediction(wide_df).summary_frame(alpha=0.05)
    wide_df['ci_var2_vs_ctrl2_lower'] = pred2['mean_ci_lower']
    wide_df['ci_var2_vs_ctrl2_upper'] = pred2['mean_ci_upper']
    wide_df['within_ci_var2'] = (wide_df['var_2'] >= wide_df['ci_var2_vs_ctrl2_lower']) & \
                                (wide_df['var_2'] <= wide_df['ci_var2_vs_ctrl2_upper'])

    # Regression 3: median of ctrls vs median of vars
    model3 = smf.ols('var_median ~ ctrl_median', data=wide_df).fit()
    pred3 = model3.get_prediction(wide_df).summary_frame(alpha=0.05)
    wide_df['ci_median_vs_median_lower'] = pred3['mean_ci_lower']
    wide_df['ci_median_vs_median_upper'] = pred3['mean_ci_upper']
    wide_df['within_ci_median'] = (wide_df['var_median'] >= wide_df['ci_median_vs_median_lower']) & \
                                  (wide_df['var_median'] <= wide_df['ci_median_vs_median_upper'])


    # --- 6. Finalize and Save ---
    final_cols = [
        'base_id', 'ctrl_1', 'ctrl_2', 'ctrl_delta', 'ctrl_var',
        'var_1', 'var_2', 'var_delta', 'var_var', 'ctrl_median', 'var_median',
        'delta_med', 'var_med', 'p_levene', 'p_ttest', 'fdr', 
        'lower_ci', 'upper_ci', 'is_significant_bootstrap',
        'ci_var1_vs_ctrl1_lower', 'ci_var1_vs_ctrl1_upper', 'within_ci_var1',
        'ci_var2_vs_ctrl2_lower', 'ci_var2_vs_ctrl2_upper', 'within_ci_var2',
        'ci_median_vs_median_lower', 'ci_median_vs_median_upper', 'within_ci_median'
    ]
    final_df = wide_df[[col for col in final_cols if col in wide_df.columns]]
    
    final_df.to_csv(output_file, index=False, sep='\t')
    print(f"\nAnalysis complete. Comprehensive table saved to: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run comprehensive analysis with classical, bootstrap, and regression tests."
    )
    parser.add_argument(
        "-i", "--input",
        type=str,
        required=True,
        help="Path to the raw replicate data file."
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        required=True,
        help="Path for the comprehensive output CSV file."
    )
    parser.add_argument(
        "--sep",
        type=str,
        default=",",
        help="Separator for the input file (Default: ',')."
    )

    args = parser.parse_args()
    analyze_data_comprehensive(args.input, args.output, args.sep)
