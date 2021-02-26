# Test the effect of different lam in contractive_loss

import tensorflow as tf
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

import keras
from keras.layers import Input, Dense
from keras.models import Model
import keras.backend as K

import numpy as np
import pickle
import datetime
import time
import sys

from sklearn import preprocessing

print('cae_lam_test starting: {}'.format(datetime.datetime.now()), file=sys.stderr)
start = time.time()

early_stopping_monitor = keras.callbacks.EarlyStopping(monitor='loss', patience=20, verbose=1)

pklfile = '../../data/all_bands.pkl'
with open(pklfile, 'rb') as fi:
    data = pickle.load(fi)

data = preprocessing.MinMaxScaler().fit_transform(np.abs(data))


def contractiveAutoencoder(feat_mat, hidden_size=800, lam=1e-5, batch_size=128, epochs=2000):

    input_size = feat_mat.shape[1]
    x = Input(shape=(input_size,))
    h = Dense(hidden_size, activation='relu', name='encoded')(x)
    r = Dense(input_size, activation='sigmoid')(h)

    ae= Model(inputs=x, outputs=r)

    def contractive_loss(y_pred, y_true):
        mse = K.mean(K.square(y_true - y_pred), axis=1)

        W = K.variable(value=ae.get_layer('encoded').get_weights()[0])  # N x N_hidden
        W = K.transpose(W)  # N_hidden x N
        h = ae.get_layer('encoded').output
        dh = h * (1 - h)  # N_batch x N_hidden

        # N_batch x N_hidden * N_hidden x 1 = N_batch x 1
        contractive = lam * K.sum(dh**2 * K.sum(W**2, axis=1), axis=1)

        return mse + contractive


    ae.compile(optimizer='adam', loss=contractive_loss, metrics=['accuracy'])
    history = ae.fit(feat_mat, feat_mat, batch_size=batch_size,epochs=epochs, callbacks = [early_stopping_monitor])  
    encoder = Model(x, h)
    

    return ae, encoder, history


lam_list = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
his = []

for val in lam_list:
    cae, cae_encoder, cae_his = contractiveAutoencoder(data, lam=val, epochs=1000)
    
    model_path = '../../models/model_cae_{}.mdl'.format(val)
    encoder_path = '../../models/encoder_cae_{}.mdl'.format(val)
    
    cae.save(model_path)
    cae_encoder.save(encoder_path)

    his.append(cae_his.history)

pklfile = '../../data/history_cae.pkl'
with open(pklfile, 'wb') as fo:
    pickle.dump(his, fo)




elapsed_time = time.strftime("%H:%M:%S", time.gmtime(time.time() - start))
print('cae_lam_test finished: {}'.format(elapsed_time), file=sys.stderr)

