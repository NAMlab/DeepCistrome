import pandas as pd
from glob import glob
pd.options.display.width=0
ensembl_entrez = pd.read_csv(filepath_or_buffer='data/coexp/mart_ensembl_to_entrez_IDs.txt', sep='\t',
                             dtype={'NCBI gene (formerly Entrezgene) ID': str})
ensembl_entrez.columns = ['Gene ID', 'Entrez ID']
ensembl_entrez_dict = dict(zip(ensembl_entrez['Entrez ID'], ensembl_entrez['Gene ID']))
print(ensembl_entrez.head())
coexp_datasets = glob('data/coexp/Ath-r.v23-11.G19783-S22408.combat_pca.subagging.z.d/*')

data = pd.read_csv(coexp_datasets[0], sep='\t', header=None, dtype={0:str})
data.columns = ['Entrez ID', coexp_datasets[0].split('/')[-1]]

gene_list = [ensembl_entrez_dict[x.split('/')[-1]] for x in coexp_datasets if x.split('/')[-1] in ensembl_entrez_dict.keys()]

for idx, f in enumerate(coexp_datasets[1:]):
    if f.split('/')[-1] in ensembl_entrez_dict.keys():
        df = pd.read_csv(f, sep='\t', header=None, dtype={0:str})
        df.columns = ['Entrez ID', f.split('/')[-1]]
        data = data.merge(df, on='Entrez ID', how='outer')

data.set_index('Entrez ID', inplace=True)
data = data.loc[data.columns.tolist(), :]
data.columns = [ensembl_entrez_dict[i] for i in data.columns]
data.index = [ensembl_entrez_dict[i] for i in data.index]
print(data.head())
data.to_csv('data/coexpression_table.tsv', sep='\t')

