#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from numpy.lib import recfunctions as rfn
from sklearn import preprocessing
from sklearn.decomposition import KernelPCA
from sklearn.cluster import SpectralClustering
import sksurv.compare


class Cluster(object):

    def __init__(self, gene, isoform, methyl, clinical):
        self.gene = gene
        self.isoform = isoform
        self.methyl = methyl
        self.clinical = clinical


    def gaussian_kernel(self, E, SK_PARA):
        m, n = np.array(E).shape
        K = np.zeros((n, n))
        for i in range(0, n):
            for j in range(0, n):
                K[i][j] = np.exp(-np.sum(np.square(np.mat(E[0:m, i]) - np.mat(E[0:m, j]))) / (2 * (SK_PARA ** 2)))
        return K


    def construct_of_SK(self, SK_PARA):
        gene = pd.read_csv(self.gene).iloc[:, 1:]
        isoform = pd.read_csv(self.isoform).iloc[:, 1:]
        methyl = pd.read_csv(self.methyl).iloc[:, 1:]

        kpca = KernelPCA(kernel="poly", degree=3)
        kpca_gene = kpca.fit_transform(preprocessing.MinMaxScaler().fit_transform(gene.T))
        kpca_isoform = kpca.fit_transform(preprocessing.MinMaxScaler().fit_transform(isoform.T))
        kpca_methyl = kpca.fit_transform(preprocessing.MinMaxScaler().fit_transform(methyl.T))

        K_gene = self.gaussian_kernel(kpca_gene.T, SK_PARA)
        K_isoform = self.gaussian_kernel(kpca_isoform.T, SK_PARA)
        K_methyl = self.gaussian_kernel(kpca_methyl.T, SK_PARA)

        return (K_gene + K_isoform + K_methyl) / 3


    def clustering(self, min_number_of_clusters, max_number_of_clusters, SK):

        clinical = np.array(
            pd.read_csv(self.clinical, usecols=['vital_status', 'days_to_death']))
        clinical_struct = rfn.unstructured_to_structured(clinical,
                                                         np.dtype([('Status', '?'), ('Survival_in_days', '<f8')]))

        result = []
        for i in range(min_number_of_clusters, max_number_of_clusters + 1):
            clustering_i = SpectralClustering(n_clusters=i, affinity='precomputed', assign_labels='discretize',
                                              random_state=0).fit(SK)
            labels_i = clustering_i.labels_
            result.append(labels_i.tolist())

        result_array = np.array(result)
        result_array_t = result_array.transpose()

        # P-value
        n_features = result_array_t.shape[1]
        for j in range(n_features):
            Xj = result_array_t[:, j:j + 1]
            res = sksurv.compare.compare_survival(clinical_struct, Xj)
            print("[cluster: %d, P-value: %.8f]" % (j + 2, res[1]))


if __name__ == '__main__':
    gene = "./datadets/LUNG/gene.csv"
    isoform = "./datadets/LUNG/isoform.csv"
    methyl = "./datadets/LUNG/methyl.csv"
    clinical = "./datadets/LUNG/clinical.csv"
    cluster = Cluster(gene, isoform, methyl, clinical)
    SK = cluster.construct_of_SK(0.64)
    cluster.clustering(2, 10, SK)
    # [cluster: 2, P - value: 0.92260736]
    # [cluster: 3, P - value: 0.08198263]
    # [cluster: 4, P - value: 0.00906062]
    # [cluster: 5, P - value: 0.05835960]
    # [cluster: 6, P - value: 0.05484691]
    # [cluster: 7, P - value: 0.16474903]
    # [cluster: 8, P - value: 0.32122098]
    # [cluster: 9, P - value: 0.16902343]
    # [cluster: 10, P - value: 0.18959946]
