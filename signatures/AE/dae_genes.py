# Autoencoder models using gene features

import tensorflow as tf
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

import keras
from keras.layers import Input, Dense
from keras.models import Model

import numpy as np
import pickle
import datetime
import time
import sys

from sklearn import preprocessing

print('dae_genes starting: {}'.format(datetime.datetime.now()), file=sys.stderr)
start = time.time()

early_stopping_monitor = keras.callbacks.EarlyStopping(monitor='loss', patience=20, verbose=1)

pklfile = '../../data/selected_gene_mat.pkl'
with open(pklfile, 'rb') as fi:
    data = pickle.load(fi)
    data = np.array(data)

data = preprocessing.MinMaxScaler().fit_transform(np.abs(data))


def denoisingAutoencoder(feat_mat, noise_factor=0.05, hidden_size=800, batch_size=128, epochs=2000):
    
    input_size = feat_mat.shape[1]
    feat_noisy = feat_mat + noise_factor * np.random.normal(loc=0.0, scale=1.0, size=feat_mat.shape) 
    feat_noisy = np.clip(feat_noisy, 0., 1.)

    x = Input(shape=(input_size,))
    h = Dense(hidden_size, activation='relu')(x)
    r = Dense(input_size, activation='sigmoid')(h)

    ae = Model(inputs=x, outputs=r)
    ae.compile(optimizer='adam', loss='mse',metrics=['accuracy'])

    history = ae.fit(feat_noisy, feat_noisy, batch_size=batch_size,epochs=epochs, callbacks = [early_stopping_monitor])
    encoder = Model(x,h)
    return ae, encoder, history


core_list = [1500, 3000, 3500, 4000]
his = []

for core_size in core_list:
    dae, dae_encoder, dae_his = denoisingAutoencoder(data, hidden_size=core_size)
    
    model_path = '../../models/Genes/model_dae_core_{}.mdl'.format(core_size)
    encoder_path = '../../models/Genes/encoder_dae_core_{}.mdl'.format(core_size)
    
    dae.save(model_path)
    dae_encoder.save(encoder_path)

    his.append(dae_his)

pklfile = '../../data/Genes/history_dae_core.pkl'
with open(pklfile, 'wb') as fo:
    pickle.dump(his, fo)
    
elapsed_time = time.strftime("%H:%M:%S", time.gmtime(time.time() - start))
print('dae_core_test finished: {}'.format(elapsed_time), file=sys.stderr)
