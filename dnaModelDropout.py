#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorflow as tf
from Bio import SeqIO
from keras.models import Sequential
from keras.layers import Dense, LSTM, Activation, Dropout
from keras import regularizers
import sys
# # GPU config
physical_devices = tf.config.list_physical_devices('GPU')
tf.config.experimental.set_memory_growth(physical_devices[0], enable=True)
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

plt.style.use('bmh')
pd.options.display.width = 0

global setName

ds = [
        ['Viruses/vaccinia_test.fasta', 'Vaccines/vaccinia_vaccine.fasta'],
        ['Viruses/pestivirus_test.fasta', 'Vaccines/pestivirus_vaccine.fasta'],
        ['Viruses/measles_test.fasta', 'Vaccines/measles_morbillivirus_vaccine.fasta'],
    ]

def float_list(seq):
    # take sequence and create a mirror list with floats mapped
    # to the base pairs, as well as the count index for the virus/vax name.
    # A = 0, T = .25, G = .5, C = .75, N = 1

    temp = []
    for j in enumerate(str(seq)):
        floating_pt_val = j[1].replace('A', '0').replace('T', '.25') \
            .replace('G', '.5').replace('C', '.75').replace('N', '.375').replace('M', '.125')
        temp.append(float(floating_pt_val))
    return np.array(temp)  # dtype float64


def amend_data(pair):
    virus_strains = []
    vaccine_strains = []
    count = 0
    train_seq = SeqIO.read(ds[int(pair)][1],
                           "fasta")

    for sequence in SeqIO.parse(
            ds[int(pair)][0],
            "fasta"
    ):
        virus_strains.append(float_list(sequence.seq))  # pass sequence and count index
        vaccine_strains.append(float_list(train_seq.seq))
        count += 1

    virus_strains = np.array(virus_strains, dtype=object)
    vaccine_strains = np.array(vaccine_strains)

    return virus_strains, vaccine_strains

def virus_model(index, subject):

    vir, vax = amend_data(subject)

    for x in range(len(vir)):
        vir[x] = np.reshape(vir[x], (1, len(vir[x]), 1))

    x1 = vir[int(index)]
    y1 = np.reshape(vax[0], (1, len(vax[0]), 1))

    # dna_model = Sequential([
    #     # 3D tensor w/ shape [batch, timesteps, feature]
    #     LSTM(200, return_sequences=True, kernel_regularizer=regularizers.l2(0.0001)),
    #     Dropout(.5),
    #     LSTM(100, kernel_regularizer=regularizers.l2(0.0001)),
    #     Dropout(.5),
    #     Dense(62, kernel_regularizer=regularizers.l2(0.0001)),
    #     Dropout(.5),
    #     Dense(1)
    # ])

    dna_model = tf.keras.models.load_model('dnaModelSel')
    dna_model.compile(
        loss='binary_crossentropy',
        optimizer='adam'  # rmsprop
                      )

    dna_model.fit(x1, y1, epochs=5)

    dna_model.save("dnaModelDropout")
    print("index:\t" + index)
    print("subject:\t" + subject)


# Subject: indexed vax/virus pair in ds, then the strain in that set
virus_model(index=sys.argv[1], subject=sys.argv[2])


# virus_model(1)
# virus_model(2)
# virus_model(4)
# virus_model(5)


# # virus_ragged = tf.ragged.constant(virus_strains)
# # Uneven input rows
#
# # tf.keras.layers.Embedding(hash_buckets, 16),
#
# # Turns positive integers (indexes) into dense vectors of fixed size.
# # e.g. [[4], [20]] -> [[0.25, 0.1], [0.6, -0.2]]
# # This layer can only be used as the first layer in a model.
#
# # tf.keras.layers.Input(shape=(199617, 1),
# #                       ragged=True),
# # # Instantiates a tensor w/ input

