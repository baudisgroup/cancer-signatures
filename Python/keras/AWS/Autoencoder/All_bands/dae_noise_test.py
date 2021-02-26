# Test the effect of different nosiy factor

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

print('dae_noise_test starting: {}'.format(datetime.datetime.now()), file=sys.stderr)
start = time.time()

early_stopping_monitor = keras.callbacks.EarlyStopping(monitor='loss', patience=20, verbose=1)

pklfile = '../../data/all_bands.pkl'
with open(pklfile, 'rb') as fi:
    data = pickle.load(fi)

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


noise_list = [0.005, 0.01, 0.025, 0.05, 0.075, 0.1]
his = []

for noise in noise_list:
    dae, dae_encoder, dae_his = denoisingAutoencoder(data, noise_factor=noise, epochs=1000)

    model_path = '../../models/model_dae_noise_{}.mdl'.format(noise)
    encoder_path = '../../models/encoder_dae_noise_{}.mdl'.format(noise)
    
    dae.save(model_path)
    dae_encoder.save(encoder_path)

    his.append(dae_his)

pklfile = '../../data/history_dae_noise.pkl'
with open(pklfile, 'wb') as fo:
    pickle.dump(his, fo)


elapsed_time = time.strftime("%H:%M:%S", time.gmtime(time.time() - start))
print('dae_noise_test finished: {}'.format(elapsed_time), file=sys.stderr)

