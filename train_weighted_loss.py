import numpy as np
from utils import train_model
import pandas as pd
import os
from tensorflow.keras.backend import clear_session

if not os.path.exists('results'):
    os.mkdir('results')

if not os.path.exists('saved_models'):
    os.mkdir('saved_models')

genome_path = 'data/genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa'
overlap_mat_path = 'data/overlap_matrix.bed'

overlap_mat = pd.read_csv(overlap_mat_path, sep='\t', index_col=0)
overlap_mat['chrom'] = [i.split(':')[0] for i in overlap_mat.index.tolist()]
chromosomes = ['1', '2', '3', '4', '5']

for chrom in chromosomes:
    train_overlap_mat = overlap_mat.copy()
    train_overlap_mat = train_overlap_mat[train_overlap_mat['chrom'] != chrom]
    print(f"Training on chromosomes {train_overlap_mat['chrom'].unique()}")
    train_overlap_mat.drop(columns=['chrom'], inplace=True)
    valid_overlap_mat = overlap_mat.copy()
    valid_overlap_mat = valid_overlap_mat[valid_overlap_mat['chrom'] == chrom]
    print(f"Validating on chromosome: {valid_overlap_mat['chrom'].unique()}")
    valid_overlap_mat.drop(columns=['chrom'], inplace=True)

    label_weights = train_overlap_mat.sum(axis=0).values/np.sum(train_overlap_mat.sum(axis=0).values)

    print('------------------------------ \n Training model \n ------------------------------ \n')
    clear_session()
    train_model(genome=genome_path, batch_size=512,
                train_overlap_mat=train_overlap_mat,
                valid_overlap_mat=valid_overlap_mat,
                window_size=250, chrom=chrom, label_weights=label_weights, control='model',
                weighted_class=True)