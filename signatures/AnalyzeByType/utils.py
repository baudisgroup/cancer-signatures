import time
import pickle
import os

import numpy as np
import pandas as pd

from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

import seaborn as sns
import matplotlib.pyplot as plt


class CancerType:

    matfile = '/Users/bogao/DataFiles/new landscape/data/reduced_gene_mat.pkl'
    labelfile = '/Users/bogao/DataFiles/new landscape/data/all_bands_label.pkl'
    dlabelfile = '/Users/bogao/DataFiles/new landscape/data/all_bands_disease_label.pkl'
    censusfile = '/Users/bogao/DataFiles/new landscape/files/Census_allThu Jan 16 08_07_54 2020.tsv'
    ampgenesfile = '/Users/bogao/DataFiles/Data/genome/reduced_high_amp_genes.tsv'
    delgenesfile = '/Users/bogao/DataFiles/Data/genome/reduced_high_del_genes.tsv'
    outputpath = '/Users/bogao/DataFiles/new landscape/files/new/'

    num_amp_genes = 750

    with open(matfile, 'rb') as fmat, open(labelfile, 'rb') as flab, open(dlabelfile, 'rb') as fdl:
            data = pickle.load(fmat)
            data = np.array(data)
            labels = pickle.load(flab)
            disease_labels = pickle.load(fdl)
            census_genes = pd.read_csv(censusfile, sep='\t')
            # census_genes['chr'], census_genes['range'] = census_genes['Genome Location'].str.split(pat=':').str
            # census_genes['start'], census_genes['end'] = census_genes['range'].str.split(pat='-').str
            # # census_genes['start'] = census_genes['start'].astype('int')
            # # census_genes['end'] = census_genes['end'].astype('int')
            # census_genes.drop('range', axis=1, inplace=True)


    # Old methods, seperate prepare function, trim disease code
    # def __init__(self, cancer_type):
    #     self.type = cancer_type
    #     self.subtypes = {}

    # def prepareData(self, trim=3, min_samples=10):
        
    #     # data of a cancer type
    #     self.data = CancerType.data[CancerType.labels == self.type]

    #     # labels of a cancer type
    #     self.labels = CancerType.disease_labels[CancerType.labels == self.type]
    #     self.labels = pd.Series(self.labels).str[:trim].values

    #     # count samples in lablels(subtypes) 
    #     unique_labels, unique_counts = np.unique(self.labels,return_counts=True)
    #     unique_labels = unique_labels[unique_counts>min_samples]

    #     # remove lablels that have too few samples
    #     self.data = self.data[np.isin(self.labels, unique_labels)]
    #     self.labels = self.labels[np.isin(self.labels, unique_labels)]

    def __init__(self, cancer_type, min_samples=10):
        self.type = cancer_type

        # data of a cancer type
        self.data = CancerType.data[CancerType.labels == self.type]

        # labels of a cancer type
        self.labels = CancerType.disease_labels[CancerType.labels == self.type]
        
        # count samples in lablels(subtypes) 
        unique_labels, unique_counts = np.unique(self.labels,return_counts=True)
        unique_labels = unique_labels[unique_counts>min_samples]

        # remove lablels that have too few samples
        self.data = self.data[np.isin(self.labels, unique_labels)]
        self.labels = self.labels[np.isin(self.labels, unique_labels)]

        self.subtypes = {}


    # Redefine disease code with a mapping dict. 
    def relabel(self, codemap):

        # remap disease code
        newlabels = np.copy(self.labels)
        for k, v in codemap.items(): 
            newlabels[self.labels==k] = v

        self.labels = newlabels



    def pcaTSNE(self, RS=1234, pca_n=50):
       
       # PCA 50
        time_start = time.time()

        pca = PCA(n_components = pca_n)
        pca_result = pca.fit_transform(self.data)

        print('PCA done! Time elapsed: {} seconds'.format(time.time()-time_start))
        print('Cumulative explained variation for 50 principal components: {}'.format(
                np.sum(pca.explained_variance_ratio_)))

        # TSNE
        time_start = time.time()
        tsne = TSNE(random_state=RS).fit_transform(pca_result)
        print('t-SNE done! Time elapsed: {} seconds'.format(time.time()-time_start))

        # plot
        df = pd.DataFrame(tsne, columns = ['ptsne1','ptsne2'])
        df['label'] = self.labels
        plt.figure(figsize=(15,10))
        sns.scatterplot(
            x="ptsne1", y="ptsne2",
            hue="label",
            palette=sns.color_palette("bright", len(set(self.labels))),
            data=df,
            legend="full",
            alpha=0.5
        )


    def groupData(self):
        self.sum = np.sum(self.data, axis=0)
        self.sum_scaled = preprocessing.minmax_scale(np.abs(self.sum))

        self.amp_genes =  pd.read_csv(CancerType.ampgenesfile, sep='\t')
        self.amp_genes = self.amp_genes.assign(cnv=self.sum[:CancerType.num_amp_genes])
        self.amp_genes = self.amp_genes.assign(cnv_scaled= self.sum_scaled[:CancerType.num_amp_genes])
        
        self.del_genes = pd.read_csv(CancerType.delgenesfile, sep='\t')
        self.del_genes = self.del_genes.assign(cnv=self.sum[CancerType.num_amp_genes:])
        self.del_genes = self.del_genes.assign(cnv_scaled= -self.sum_scaled[CancerType.num_amp_genes:])

        f, axes = plt.subplots(1, 2, figsize=(10,5))

        sns.distplot(self.amp_genes['cnv'], label='amp', ax=axes[0])
        sns.distplot(self.del_genes['cnv'], label='del', ax=axes[1])

    # def censusOverlap(self, target):
    #     return CancerType.census_genes[CancerType.census_genes['Gene Symbol'].isin(target['symbol'].values)]

    def censusOverlap(self, target):
        return target[target['symbol'].isin(CancerType.census_genes['Gene Symbol'].values)]
        # return CancerType.census_genes[CancerType.census_genes['Gene Symbol'].isin(target['symbol'].values)]    

    def significantGenes(self, amp_genes, del_genes, thresh_values):
        amp_genes = amp_genes[amp_genes['cnv']>thresh_values[1]].sort_values('band')
        del_genes = del_genes[del_genes['cnv']<thresh_values[3]].sort_values('band')
        genes = amp_genes['symbol'].dropna().tolist() + del_genes['symbol'].dropna().tolist()

        high_amp_genes = amp_genes[amp_genes['cnv']>=thresh_values[0]]
        low_amp_genes = amp_genes[amp_genes['cnv']<thresh_values[0]]

        high_del_genes = del_genes[del_genes['cnv']<=thresh_values[2]]
        low_del_genes = del_genes[del_genes['cnv']>thresh_values[2]]

        high_amp_census = self.censusOverlap(high_amp_genes)
        low_amp_census = self.censusOverlap(low_amp_genes)
        amp_census = self.censusOverlap(amp_genes)

        high_del_census = self.censusOverlap(high_del_genes)
        low_del_census = self.censusOverlap(low_del_genes)
        del_census = self.censusOverlap(del_genes)      

        census = amp_census['symbol'].dropna().tolist() + del_census['symbol'].dropna().tolist()

        return {'amp_genes':amp_genes, 'del_genes':del_genes, 'genes':genes, 
                'high_amp_genes':high_amp_genes, 'low_amp_genes':low_amp_genes, 
                'high_del_genes':high_del_genes, 'low_del_genes':low_del_genes,
                'high_amp_census':high_amp_census, 'low_amp_census':low_amp_census, 'amp_census':amp_census,
                'high_del_census':high_del_census, 'low_del_census':low_del_census, 'del_census':del_census,
                'census':census}



    # thresh_values = [amp_high, amp_low, del_high, del_low]
    def analyze(self, thresh_values):
        # self.amp_genes = self.amp_genes[self.amp_genes['cnv']>thresh_values[1]].sort_values('band')
        # self.del_genes = self.del_genes[self.del_genes['cnv']<thresh_values[3]].sort_values('band')
        # self.genes = self.amp_genes['symbol'].dropna().tolist() + self.del_genes['symbol'].dropna().tolist()

        # self.high_amp_genes = self.amp_genes[self.amp_genes['cnv']>=thresh_values[0]]
        # self.low_amp_genes = self.amp_genes[self.amp_genes['cnv']<thresh_values[0]]

        # self.high_del_genes = self.del_genes[self.del_genes['cnv']<=thresh_values[2]]
        # self.low_del_genes = self.del_genes[self.del_genes['cnv']>thresh_values[2]]

        # self.high_amp_census = self.censusOverlap(self.high_amp_genes)
        # self.low_amp_census = self.censusOverlap(self.low_amp_genes)
        # self.amp_census = self.censusOverlap(self.amp_genes)

        # self.high_del_census = self.censusOverlap(self.high_del_genes)
        # self.low_del_census = self.censusOverlap(self.low_del_genes)
        # self.del_census = self.censusOverlap(self.del_genes)

        # self.census = self.amp_census['Gene Symbol'].dropna().tolist() + self.del_census['Gene Symbol'].dropna().tolist()

        res = self.significantGenes(self.amp_genes, self.del_genes, thresh_values)
        self.amp_genes = res['amp_genes']
        self.del_genes = res['del_genes']
        self.genes = res['genes']

        self.high_amp_genes = res['high_amp_genes']
        self.low_amp_genes = res['low_amp_genes']

        self.high_del_genes = res['high_del_genes']
        self.low_del_genes = res['low_del_genes']

        self.high_amp_census = res['high_amp_census']
        self.low_amp_census = res['low_amp_census']
        self.amp_census = res['amp_census']

        self.high_del_census = res['high_del_census']
        self.low_del_census = res['low_del_census']
        self.del_census = res['del_census']

        self.census = res['census']


    def countData(self):
        print('{}\t{}'.format('amp_genes', self.amp_genes.shape))
        print('{}\t{}'.format('del_genes', self.del_genes.shape))
        print('{}\t{}'.format('genes', len(self.amp_genes)))

        print('{}\t{}'.format('high_amp_genes', self.high_amp_genes.shape))
        print('{}\t{}'.format('low_amp_genes', self.low_amp_genes.shape))

        print('{}\t{}'.format('high_del_genes', self.high_del_genes.shape))
        print('{}\t{}'.format('low_del_genes', self.low_del_genes.shape))

        print('{}\t{}'.format('high_amp_census', self.high_amp_census.shape))
        print('{}\t{}'.format('low_amp_census', self.low_amp_census.shape))
        print('{}\t{}'.format('amp_census', self.amp_census.shape))

        print('{}\t{}'.format('high_del_census', self.high_del_census.shape))
        print('{}\t{}'.format('low_del_census', self.low_del_census.shape))
        print('{}\t{}'.format('del_census', self.del_census.shape))

        print('{}\t{}'.format('census', len(self.census)))




    def countSubtypes(self):
        print(np.array(np.unique(self.labels, return_counts=True)).T)

    def dumpFiles(self):

        path = CancerType.outputpath + self.type + '/'
        os.makedirs(path, exist_ok=True)
        path = path + self.type.lower()

        genefile = path + '_genes.txt'
        with open(genefile, 'w') as fo:
            for i in self.genes:
                print(i, file=fo)

        censusfile = path + '_census.txt'
        with open(censusfile, 'w') as fo:
            for i in self.census:
                print(i, file=fo)

        highgenefile = path + '_genes_high.txt'
        with open(highgenefile, 'w') as fo:
            for i in self.high_amp_genes['symbol'].dropna().values:
                print(i, file=fo)
            for i in self.high_del_genes['symbol'].dropna().values:
                print(i, file=fo)

        lowgenefile = path + '_genes_low.txt'
        with open(lowgenefile, 'w') as fo:
            for i in self.low_amp_genes['symbol'].dropna().values:
                print(i, file=fo)
            for i in self.low_del_genes['symbol'].dropna().values:
                print(i, file=fo)

        highcensusfile = path + '_census_high.txt'
        with open(highcensusfile, 'w') as fo:
            for i in self.high_amp_census['symbol'].dropna().values:
                print(i, file=fo)
            for i in self.high_del_census['symbol'].dropna().values:
                print(i, file=fo)

        lowcensusfile = path + '_census_low.txt'
        with open(lowcensusfile, 'w') as fo:
            for i in self.low_amp_census['symbol'].dropna().values:
                print(i, file=fo)
            for i in self.low_del_census['symbol'].dropna().values:
                print(i, file=fo)


        self.high_amp_genes.to_csv(path + '_high_amp_genes.tsv', sep='\t', index=False)
        self.low_amp_genes.to_csv(path + '_low_amp_genes.tsv', sep='\t', index=False)

        self.high_del_genes.to_csv(path + '_high_del_genes.tsv', sep='\t', index=False)
        self.low_del_genes.to_csv(path + '_low_del_genes.tsv', sep='\t', index=False)

        self.high_amp_census.to_csv(path + '_high_amp_census.tsv', sep='\t', index=False)
        self.low_amp_census.to_csv(path + '_low_amp_census.tsv', sep='\t', index=False)

        self.high_del_census.to_csv(path + '_high_del_census.tsv', sep='\t', index=False)
        self.low_del_census.to_csv(path + '_low_del_census.tsv', sep='\t', index=False)

        self.amp_genes.to_csv(path + '_amp_genes.tsv', sep='\t', index=False)
        self.del_genes.to_csv(path + '_del_genes.tsv', sep='\t', index=False)

        self.amp_census.to_csv(path + '_amp_census.tsv', sep='\t', index=False)
        self.del_census.to_csv(path + '_del_census.tsv', sep='\t', index=False)





    def prepareSubtype(self, subtype):
        subdata = self.data[self.labels == subtype]
        subsum = np.sum(subdata, axis=0)
        subsum_scaled = preprocessing.minmax_scale(np.abs(subsum))

        amp_genes = pd.read_csv(CancerType.ampgenesfile, sep='\t')
        amp_genes = amp_genes.assign(cnv = subsum[:CancerType.num_amp_genes])
        amp_genes = amp_genes.assign(cnv_scaled= subsum_scaled[:CancerType.num_amp_genes])

        del_genes = pd.read_csv(CancerType.delgenesfile, sep='\t')
        del_genes = del_genes.assign(cnv = subsum[CancerType.num_amp_genes:])
        del_genes = del_genes.assign(cnv_scaled= -subsum_scaled[CancerType.num_amp_genes:])

        f, axes = plt.subplots(1, 2, figsize=(10,5))

        sns.distplot(amp_genes['cnv'], label='amp', ax=axes[0])
        sns.distplot(del_genes['cnv'], label='del', ax=axes[1])

        self.subtypes[subtype] = {'amp_genes':amp_genes, 'del_genes':del_genes}


    def analyzeSubtype(self, subtype, thresh_values):

        self.subtypes[subtype] = self.significantGenes( self.subtypes[subtype]['amp_genes'], 
                                                        self.subtypes[subtype]['del_genes'], 
                                                        thresh_values)
    
    def countSubtypeData(self, subtype):
        for k, i in self.subtypes[subtype].items():
            print('{}\t{}'.format(k, len(i)))


    def dumpSubtypeFiles(self, subtype):
        for k, i in self.subtypes[subtype].items():
            path = '{}{}/{}'.format(CancerType.outputpath, self.type, subtype)
            os.makedirs(path, exist_ok=True)
            path = '{}/{}_{}_'.format(path, self.type.lower(), subtype)

            if type(i) is list:
                with open(path + k + '.tsv', 'w') as fo:
                    for line in i:
                        print(line, file=fo)
            else:
                i.to_csv(path + k + '.tsv', sep='\t', index=False)

    def dumpSubtypeCounts(self):
        counts = pd.DataFrame(np.array(np.unique(self.labels, return_counts=True)).T, columns=['name','count'])
        path = '{}{}/{}_counts.tsv'.format(CancerType.outputpath, self.type, self.type.lower())
        counts.to_csv(path, index=False, sep='\t')
















