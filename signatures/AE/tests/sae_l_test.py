# Test the effect of different L in regularizer

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

print('sae_l_test starting: {}'.format(datetime.datetime.now()), file=sys.stderr)
start = time.time()

early_stopping_monitor = keras.callbacks.EarlyStopping(monitor='loss', patience=20, verbose=1)

pklfile = '../../data/all_bands.pkl'
with open(pklfile, 'rb') as fi:
    data = pickle.load(fi)

data = preprocessing.MinMaxScaler().fit_transform(np.abs(data))


def sparseAutoencoder(feat_mat, l=1e-5, hidden_size=800, batch_size=128, epochs=2000):
    
    input_size = feat_mat.shape[1]


    x = Input(shape=(input_size,))
    h = Dense(hidden_size, activation='relu', activity_regularizer=keras.regularizers.l1(l))(x)
    r = Dense(input_size, activation='sigmoid')(h)

    ae = Model(inputs=x, outputs=r)
    ae.compile(optimizer='adam', loss='mse',metrics=['accuracy'])

    history = ae.fit(feat_mat, feat_mat, batch_size=batch_size,epochs=epochs, callbacks = [early_stopping_monitor])
    encoder = Model(x,h)
    return ae, encoder, history


l_list = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
his = []

for val in l_list:
    sae, sae_encoder, sae_his = sparseAutoencoder(data, l=val, epochs=1000)
    
    model_path = '../../models/model_sae_{}.mdl'.format(val)
    encoder_path = '../../models/encoder_sae_{}.mdl'.format(val)
    
    sae.save(model_path)
    sae_encoder.save(encoder_path)

    his.append(sae_his)

pklfile = '../../data/history_sae.pkl'
with open(pklfile, 'wb') as fo:
    pickle.dump(his, fo)



elapsed_time = time.strftime("%H:%M:%S", time.gmtime(time.time() - start))
print('sae_l_test finished: {}'.format(elapsed_time), file=sys.stderr)