import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import matplotlib.pyplot as plt
np.random.seed(42)
sns.set_context(context="paper", rc={"font.size":14,"axes.titlesize":14,"axes.labelsize":14}, font_scale=1.4)

# Load coexpression matrix and regulatory cluster information
coexp_z = pd.read_csv(filepath_or_buffer='data/coexpression_table.tsv', sep='\t', index_col=0)
coexp_z = coexp_z[~coexp_z.index.duplicated(keep='first')]
clusters = pd.read_csv(filepath_or_buffer="data/prom_term_predictions.csv")
clusters.rename(mapper={'Unnamed: 0': 'gene_id'}, axis='columns', inplace=True)
clusters = clusters[['gene_id', 'cluster']]
clusters = clusters[clusters['cluster'] != 0]

# keep only shared genes
shared_genes = list(set(coexp_z.index).intersection(clusters['gene_id'].unique()))
coexp_z = coexp_z.loc[shared_genes, shared_genes]
clusters = clusters.set_index('gene_id').loc[shared_genes]
clusters['gene_id'] = clusters.index.tolist()
gene_to_idx = {gene: i for i, gene in enumerate(clusters['gene_id'].tolist())}
clusters = clusters.groupby('cluster')['gene_id'].apply(list).to_dict()
coexp_z = coexp_z.values

def mean_coexp(cluster_genes):
    idx = [gene_to_idx[g] for g in cluster_genes]
    # get upper triangle indices (i<j) for pairs
    tri_upper = np.triu_indices(len(idx), k=1)
    # extract the z-scores for these pairs
    pair_z = coexp_z[np.ix_(idx, idx)][tri_upper]
    return pair_z.mean()


# Permutation test
def cluster_pval(cluster_genes, n_perm=1000):
    k = len(cluster_genes)
    if k < 2:
        return np.nan, np.nan, np.nan
    obs_mean = mean_coexp(cluster_genes)

    # null distribution
    null_means = []
    all_genes = list(gene_to_idx.keys())
    for _ in range(n_perm):
        sample_genes = np.random.choice(all_genes, size=k, replace=False)
        null_means.append(mean_coexp(sample_genes))

    null_means = np.array(null_means)
    # p-value: fraction of null >= observed
    pval = (np.sum(null_means >= obs_mean) + 1) / (n_perm + 1)

    # signed standardized effect size
    null_mean = null_means.mean()
    null_sd = null_means.std(ddof=1)

    effect_size = (obs_mean - null_mean) / null_sd
    return pval, obs_mean, effect_size


# Compute p-values for all clusters
cluster_ids = []
mean_vals = []
pvals = []
effect_sizes = []
for cl_id, genes in clusters.items():
    print(cl_id)
    cluster_ids.append(f'c{cl_id}')
    p, obs_mean, effect_s = cluster_pval(genes, n_perm=1000)
    mean_vals.append(obs_mean)
    pvals.append(p)
    effect_sizes.append(effect_s)

# Multiple testing correction
reject, pvals_corrected, _, _ = multipletests(pvals, method='fdr_bh')

# Output results including mean coexpression
results = pd.DataFrame({
    'cluster': cluster_ids,
    'mean coexpression logit score': mean_vals,
    'pval': pvals,
    'pval_adj': pvals_corrected,
    'significant': reject,
    'effect size': effect_sizes
})

# Plotting---------------------
# Add -log10(FDR-adjusted p-value)
results['-log10(adj-Pvalue)'] = -np.log10(results['pval_adj'])

results.to_csv(path_or_buf="results/cluster_coexpression_significance.csv", index=False, sep='\t')
print(results.head(15))

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 4))
sns.scatterplot(
    data=results,
    x='mean coexpression logit score',
    y='-log10(adj-Pvalue)',
    hue='significant',
    palette={True:'#EB4C4C', False:'#BFC9D1'},
    s=120,
    ax=ax
)

# Optional: horizontal line for FDR = 0.05
ax.axhline(-np.log10(0.05), color='#7EACB5', linestyle='--', label='FDR 0.05')
ax.yaxis.grid(True, alpha=0.3)
ax.xaxis.grid(True, alpha=0.3)
ax.set_axisbelow(True)
ax.spines[['right', 'top']].set_visible(False)
fig.tight_layout()

plt.legend(title='Significant')
plt.savefig(f"results/Figures/Coexpression_cluster_analysis.svg", bbox_inches='tight',
                dpi=300, format='svg')
plt.show()