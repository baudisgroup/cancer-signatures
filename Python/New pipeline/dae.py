# Test the effect of different core size 

# import tensorflow as tf
# tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

import os
os.environ["KERAS_BACKEND"] = "plaidml.keras.backend"

import keras
from keras.layers import Input, Dense
from keras.models import Model

import numpy as np
import pickle
import datetime
import time
import sys

from sklearn import preprocessing



print('dae_core_test starting: {}'.format(datetime.datetime.now()), file=sys.stderr)
start = time.time()

early_stopping_monitor = keras.callbacks.EarlyStopping(monitor='loss', patience=20, verbose=1)

pklfile = '/Users/bogao/DataFiles/new landscape/data/band2gene_selected_gene_mat_array.pkl'
# pklfile = '/Users/bogao/DataFiles/new landscape/data/selected_gene_mat.pkl'
with open(pklfile, 'rb') as fi:
    data = pickle.load(fi)
    # data = np.array(data)

data = preprocessing.MinMaxScaler().fit_transform(np.abs(data))

# pklfile = '/Users/bogao/DataFiles/new landscape/data/all_bands_label.pkl'
# with open(pklfile, 'rb') as fi:
#     labels = pickle.load(fi)
# data = data[labels != 'Others']

def denoisingAutoencoder(feat_mat, noise_factor=0.05, hidden_size=800, batch_size=128, epochs=1000):
    
    input_size = feat_mat.shape[1]
    feat_noisy = feat_mat + noise_factor * np.random.normal(loc=0.0, scale=1.0, size=feat_mat.shape) 
    feat_noisy = np.clip(feat_noisy, 0., 1.)

    x = Input(shape=(input_size,))

    # hidden_encoder_1 = Dense(int(input_size/4), activation='relu')(x)
    # hidden_encoder_2 = Dense(int(input_size/4), activation='relu')(hidden_encoder_1)
    # hidden_encoder_3 = Dense(int(input_size/4), activation='relu')(hidden_encoder_2)

    h = Dense(hidden_size, activation='relu')(x)

    # hidden_decoder_1 = Dense(int(input_size/4), activation='relu')(h)
    # hidden_decoder_2 = Dense(int(input_size/4), activation='relu')(hidden_decoder_1)
    # hidden_decoder_3 = Dense(int(input_size/4), activation='relu')(hidden_decoder_2)

    r = Dense(input_size, activation='sigmoid')(h)

    ae = Model(inputs=x, outputs=r)
    ae.compile(optimizer='adam', loss='mse',metrics=['accuracy'])

    history = ae.fit(feat_noisy, feat_noisy, batch_size=batch_size,epochs=epochs, callbacks = [early_stopping_monitor])
    encoder = Model(x,h)
    return ae, encoder, history


core_list = [2000]

for core_size in core_list:
    dae, dae_encoder, dae_his = denoisingAutoencoder(data[:,:2000], hidden_size=core_size)
    
    model_path = '/Users/bogao/DataFiles/new landscape/models/new_pipeline_dae_genes_model_{}.mdl'.format(core_size)
    encoder_path = '/Users/bogao/DataFiles/new landscape/models/new_pipeline_dae_genes_encoder_{}.mdl'.format(core_size)
    his_path = '/Users/bogao/DataFiles/new landscape/models/new_pipline_dae_history_{}.pkl'.format(core_size)
    
    dae.save(model_path)
    dae_encoder.save(encoder_path)

    with open(his_path, 'wb') as fo:
        pickle.dump(dae_his, fo)
    
elapsed_time = time.strftime("%H:%M:%S", time.gmtime(time.time() - start))
print('dae_core_test finished: {}'.format(elapsed_time), file=sys.stderr)
