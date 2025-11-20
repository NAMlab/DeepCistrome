from pybedtools import BedTool
import os
from tensorflow.keras.utils import Sequence
import math
from pyfaidx import Fasta
import numpy as np
from typing import Any
import tensorflow as tf
import tensorflow.keras.layers as kl
from tensorflow.keras import Model, optimizers
from tensorflow.keras.metrics import AUC
from deeplift.dinuc_shuffle import dinuc_shuffle
import pandas as pd
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint, CSVLogger
print(f"Number of GPUs: { len(tf.config.list_physical_devices('GPU'))}")


def one_hot_encode(sequence: str,
                   alphabet: str = 'ACGT',
                   neutral_alphabet: str = 'N',
                   neutral_value: Any = 0,
                   dtype=np.float32) -> np.ndarray:
    """One-hot encode sequence."""
    def to_uint8(string):
        return np.frombuffer(string.encode('ascii'), dtype=np.uint8)
    hash_table = np.zeros((np.iinfo(np.uint8).max, len(alphabet)), dtype=dtype)
    hash_table[to_uint8(alphabet)] = np.eye(len(alphabet), dtype=dtype)
    hash_table[to_uint8(neutral_alphabet)] = neutral_value
    hash_table = hash_table.astype(dtype)
    return hash_table[to_uint8(sequence)]


def create_genome_bins(genome, bin_size):
    os.system(f"samtools faidx {genome}")
    genome_bins = BedTool().window_maker(g=f"{genome}.fai",
                                         w=bin_size, s=bin_size).filter(lambda r: r.chrom in ['1', '2', '3', '4', '5'])
    genome_bins.saveas(f"{genome.split('_')[0]}_genome_bins.bed")


def peak_overlap_genome_bins(genome, frac, bin_size=250):
    create_genome_bins(genome=genome, bin_size=bin_size)
    genome_bins_bed = f"{genome.split('_')[0]}_genome_bins.bed"
    os.system(f"bedtools intersect -a {genome_bins_bed} -b data/peaks/*.bed -C -filenames -F {frac} > data/overlap.bed")


class InputGenerator(Sequence):
    def __init__(self, genome, bins_list, batch_size, overlap_mat, window_size, control):
        self.genome = Fasta(genome, as_raw=True, read_ahead=10000, sequence_always_upper=True)
        self.bins_list = bins_list
        self.batch_size = batch_size
        self.overlap_mat = overlap_mat
        self.window_size = window_size
        self.control = control

    def __len__(self):
        return math.ceil(len(self.bins_list) / self.batch_size)

    def __getitem__(self, idx):
        batch_x = self.bins_list[idx * self.batch_size:(idx + 1) * self.batch_size]
        x = []
        y = []

        for bin_coord in batch_x:
            chrom, start, end = bin_coord.split(':')
            seq = one_hot_encode(self.genome[chrom][int(start):int(end)])
            if self.control == 'si':
                np.random.shuffle(seq)
            elif self.control == 'di':
                seq = dinuc_shuffle(seq)

            if seq.shape[0] == self.window_size:
                x.append(seq)
                y.append(self.overlap_mat.loc[bin_coord, :].tolist())

        x, y = np.array(x), np.array(y)
        return x, y


def compiled_model(window_size, num_tf_fam, label_weights):
    inputs = kl.Input(shape=(window_size, 4))
    x = kl.Conv1D(filters=256, kernel_size=21, padding='same')(inputs)
    x = kl.BatchNormalization()(x)
    x = kl.Activation('relu')(x)
    x = kl.MaxPooling1D(2)(x)

    conv_filters, conv_kernel_sizes = [60, 60, 120], [7, 7, 5]
    for n_filters, kernel_size in zip(conv_filters, conv_kernel_sizes):
        x = kl.Conv1D(filters=n_filters, kernel_size=kernel_size, padding='same')(x)
        x = kl.BatchNormalization()(x)
        x = kl.Activation('relu')(x)
        x = kl.MaxPooling1D(2)(x)

    x = kl.Flatten()(x)
    for units in [256, 256]:
        x = kl.Dense(units=units)(x)
        x = kl.BatchNormalization()(x)
        x = kl.Activation('relu')(x)
        x = kl.Dropout(0.4)(x)

    x = kl.Dense(units=num_tf_fam)(x)
    x = kl.Activation('sigmoid')(x)
    model = Model(inputs=inputs, outputs=x)
    model.compile(loss='binary_crossentropy',
                  metrics=[AUC(multi_label=True, curve='ROC', name='auROC'),
                           AUC(multi_label=True, curve='PR', name='auPR'),
                           AUC(multi_label=True, curve='PR', name='weighted_auPR', label_weights=label_weights)
                           ],
                  optimizer=optimizers.Adam(lr=0.002))
    model.summary()
    return model


def train_model(genome, batch_size, train_overlap_mat, valid_overlap_mat, window_size, chrom, label_weights,
                control, weighted_class=False):
    if weighted_class:
        model_save_name = f"saved_models/model_chrom_{chrom}_{control}_weighted.h5"
        cvs_log_save_name = f"results/csv_{chrom}_{control}_weighted.log"
    else:
        model_save_name = f"saved_models/model_chrom_{chrom}_{control}.h5"
        cvs_log_save_name = f"results/csv_{chrom}_{control}.log"

    if not os.path.exists(model_save_name) and not os.path.exists(cvs_log_save_name):
        train_gen = InputGenerator(genome=genome, bins_list=train_overlap_mat.index.tolist(), batch_size=batch_size,
                                   overlap_mat=train_overlap_mat, window_size=window_size, control=control)
        valid_gen = InputGenerator(genome=genome, bins_list=valid_overlap_mat.index.tolist(), batch_size=batch_size,
                                   overlap_mat=valid_overlap_mat, window_size=window_size, control=None)

        model = compiled_model(window_size=window_size, num_tf_fam=train_overlap_mat.shape[1],
                               label_weights=label_weights)

        if weighted_class:
            weights = (1 - label_weights/max(label_weights))+0.02
            class_weights = {k:v for k,v in enumerate(weights)}
            print('These are the class weights:\n', class_weights)
            model.fit(train_gen, validation_data=valid_gen, epochs=50,
                      callbacks=[
                          ModelCheckpoint(filepath=model_save_name, monitor='val_loss', save_best_only=True, verbose=1),
                          EarlyStopping(patience=5),
                          CSVLogger(cvs_log_save_name)
                      ], class_weight=class_weights)
        else:
            model.fit(train_gen, validation_data=valid_gen, epochs=50,
                      callbacks=[
                          ModelCheckpoint(filepath=model_save_name, monitor='val_loss', save_best_only=True, verbose=1),
                          EarlyStopping(patience=5),
                          CSVLogger(cvs_log_save_name)
                      ])
