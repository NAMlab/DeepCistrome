import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
pd.options.display.width=0
data_coexp = pd.read_csv(filepath_or_buffer='data/coexpression_table.tsv', sep='\t', index_col=0)
cluster_idx_to_name = pd.read_csv(filepath_or_buffer='data/cluster_ids_cluster_name.csv')
predicted_clusters = pd.read_csv(filepath_or_buffer="data/prom_term_predictions.csv")
predicted_clusters.rename({'Unnamed: 0': 'gene_id'}, axis='columns', inplace=True)
predicted_clusters = predicted_clusters[['gene_id', 'cluster', 'cluster_col']]
coexp_genes = data_coexp.index.tolist()
print(data_coexp.shape)
print(data_coexp.head())
print(cluster_idx_to_name.head())
print(predicted_clusters.head())
cluster_to_cluster_id_dict = {k:v for k,v in cluster_idx_to_name.values}
final_result, final_result_adjusted = [], []
num_iterations = 100
for cluster_name, cluster_genes in predicted_clusters.groupby('cluster'):
    print(f'Processing: {cluster_name}')
    genes = cluster_genes['gene_id'].tolist()
    genes_found = [x for x in genes if x in coexp_genes]
    logit_scores_genes_found = data_coexp.loc[genes_found, genes_found].values
    logit_scores_genes_found = list(logit_scores_genes_found[np.triu_indices_from(logit_scores_genes_found)])
    result = [cluster_name, cluster_to_cluster_id_dict[cluster_name], len(genes_found)]
    result_adjusted = [cluster_name, cluster_to_cluster_id_dict[cluster_name], len(genes_found)]
    for idx in range(1, num_iterations+1):
        logit_scores_random_genes = data_coexp.sample(n=len(genes_found), replace=False)
        logit_scores_random_genes = logit_scores_random_genes.loc[:, logit_scores_random_genes.index.tolist()].values
        logit_scores_random_genes = list(logit_scores_random_genes[np.triu_indices_from(logit_scores_random_genes)])
        mann_wu = mannwhitneyu(x=logit_scores_genes_found, y=logit_scores_random_genes, nan_policy='omit',
                               alternative='two-sided')
        result.append(mann_wu.pvalue)
        result_adjusted.append(mann_wu.pvalue * num_iterations)
    final_result.append(result)
    final_result_adjusted.append(result_adjusted)

columns = ['cluster_number', 'cluster_id', 'number of genes']
num_test = [f'test_{i}' for i in range(1, num_iterations+1)]
columns.extend(num_test)
final_result = pd.DataFrame(final_result, columns=columns)
final_result_adjusted = pd.DataFrame(final_result_adjusted, columns=columns)
print(final_result.head(15))
print(final_result_adjusted.head(15))
final_result_adjusted.to_csv(path_or_buf='results/coexpression_cluster_to_random_sig_table_adjusted.tsv', sep='\t', index=False)
final_result.to_csv(path_or_buf='results/coexpression_cluster_to_random_sig_table.tsv', sep='\t', index=False)




