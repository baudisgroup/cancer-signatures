import numpy as np
import pandas as pd

# import smote_variants as sv

from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn import utils
from sklearn.model_selection import StratifiedShuffleSplit

import matplotlib.pyplot as plt
import seaborn as sns


def cross_validation(data, target, spliter, sampler, model):
    
    result = []
    for train_index, test_index in spliter.split(data, target):
        X_train, X_test = data[train_index], data[test_index]
        y_train, y_test = target[train_index], target[test_index]

        if sampler == None:
            X_samp, y_samp = X_train, y_train
        else:
            X_samp, y_samp = sampler.sample(X_train, y_train)

        model.fit(X_samp, y_samp)
        y_pred = model.predict(X_test)

        result.append([y_test, y_pred])
    return result

def evaluate_global(results):
    for r,i in zip(results, range(len(results))):
        print('-------- {} --------'.format(i))
        print('Accuracy: {:.4f}'.format(metrics.accuracy_score(r[0], r[1])))
        print('Precision: {:.4f}'.format(metrics.precision_score(r[0], r[1], average='macro')))
        print('Recall: {:.4f}'.format(metrics.recall_score(r[0], r[1], average='macro')))
        print('F1: {:.4f}'.format(metrics.f1_score(r[0], r[1], average='macro'))) 

def evaluate_classes(result, names, sort='F1-score', ascending=False):
    prfs = metrics.precision_recall_fscore_support(result[0], result[1], average=None)
    metrics_table = pd.DataFrame({'Label':names, 'Precision': prfs[0], 'Recall':prfs[1], 'F1-score':prfs[2]})
    return metrics_table.sort_values(sort, ascending=ascending)

def plot_confusion_matrix(result, names, scale=False):
    cm = metrics.confusion_matrix(result[0], result[1])
    if scale == True:
        cm = cm / cm.sum(axis=1, keepdims=True)
    
    plt.figure(figsize=(11,10))
    ax = sns.heatmap(cm , cmap="YlGnBu", xticklabels=names, yticklabels=names)
    blim, tlim = ax.get_ylim()
    ax = ax.set_ylim(blim+0.5, tlim-0.5)

def false_predictions(result, names, label, sort='sum', ascending=False):
    code = np.where(names == label)[0][0]
    trues = result[0]
    preds = result[1]
    
    # fasle negative
    wrongs = []
    for t, p in zip(trues, preds):
        if t == code and p != code:
            wrongs.append(p)
    fn = pd.DataFrame(np.unique([names[x] for x in wrongs], return_counts=True)).transpose().rename(columns={0:'name', 1:'false_neg'})
    
    # false postive
    wrongs = []
    for t, p in zip(trues, preds):
        if t != code and p == code:
            wrongs.append(t)
    fp = pd.DataFrame(np.unique([names[x]  for x in wrongs], return_counts=True)).transpose().rename(columns={0:'name', 1:'false_pos'})
    
    # combine
    false_preds = pd.merge(fn, fp, how='outer', on='name')
    false_preds = false_preds.fillna(0)
    false_preds['sum'] = false_preds['false_neg'] + false_preds['false_pos']
    
    false_preds = false_preds.sort_values(sort, ascending=ascending)
    return false_preds

def organ_labels(result, names):
    trues = result[0]
    preds = result[1]
    
    label2organ = {}
    for v,i in zip(names, range(len(names))):
        label2organ[i] = v.split(' ')[0]
        
    trues_org = [label2organ[x] for x in trues]
    preds_org = [label2organ[x] for x in preds]
    return [trues_org, preds_org]

def get_organs(names):
    organs = []
    for i in names:
        organs.append(i.split(' ')[0])
    return np.unique(organs)

def under_sample(data, label, samples):
    
    target_idx = np.where(data['label'] == label)[0]

    samples = target_idx.shape[0] - samples
    exclude_idx = utils.resample(target_idx, n_samples=samples, replace=False)

    mask = np.ones(data['data'].shape[0], bool)
    mask[exclude_idx] = False

    us_data = data['data'][mask]
    us_target = data['target'][mask]
    
    return [us_data, us_target]