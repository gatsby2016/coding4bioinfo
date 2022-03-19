import sys
import os
from matplotlib.pyplot import axis

import math
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn import preprocessing
from sklearn.metrics.pairwise import pairwise_distances

class SC3(object):
    def __init__(self, data_root=None, anno_root=None):
        print(">>>SC3 init")
        assert data_root is not None, "please pass data root!"
        self.adata = sc.read_csv(data_root).transpose() # csv文件中行为基因列为cell sc读进来行列倒置了因此需要转置
        # self.adata.X # matrix val
        # self.adata.obs # cell name, panda dataframe
        # self.adata.var # gene name, panda dataframe

        if anno_root is not None:
            self.anno = sc.read(anno_root) 
    
    def sc3_one_go_run(self):
        print(">>>SC3 one-go running")
        self.sc3_preprocess()
        est_k = self.sc3_estimate_k()
        self.sc3_calc_dists()
        pass

    # self.adata has been gene filter and log2(X+1) trans.
    def sc3_preprocess(self, pct_dropout_min=0.1, pct_dropout_max=0.9):
        print(">>>SC3 preprocessing")
        num_cell = self.adata.shape[0]

        sc.pp.filter_genes(self.adata, min_cells=pct_dropout_min*num_cell+1) # align to R code, not equal 
        sc.pp.filter_genes(self.adata, max_cells=pct_dropout_max*num_cell-1) # gene filter between (min, max)
        sc.pp.log1p(self.adata, base=2) # log norm
        print(">>>Matrix shape:", self.adata.shape)
        
    # Estimate the optimal k for k-means clustering: align to R code
    def sc3_estimate_k(self):
        print(">>>SC3 estimate k val")
        num_cell, num_gene = self.adata.shape
        muTW = pow((math.sqrt(num_gene - 1) + math.sqrt(num_cell)), 2)
        sigmaTW = (math.sqrt(num_gene - 1) + math.sqrt(num_cell)) * pow((1/math.sqrt(num_gene - 1) + 1/math.sqrt(num_cell)), 1/3)
        bd = 3.273 * sigmaTW + muTW  # 3.2730 is the p=0.001 percentile point for the Tracy-Widom distribution

        x = preprocessing.scale(self.adata.X, axis=1)
        sigmaHatNaive = np.matmul(x, x.transpose())
        
        # compute eigenvalues and return the amount which falls above the bound
        evals, _ = np.linalg.eig(sigmaHatNaive)
        # k = 0
        # for i in range(len(evals)):
        #     if evals[i] > bd:
        #         k += 1
        k = (evals > bd).sum()
        return k

    # Calculate a distance matrix
    #  Distance between the cells are calculated using `Euclidean, Pearson and Spearman` metrics to construct.
    def sc3_calc_dists(self, Euclidean=True, Pearson=True, Spearman=True):
        print(">>>SC3 cal dists")
        data_matrix = self.adata.X

        if Euclidean:
            print(">>>SC3 cal Euclidean dists")
            pairwise_distances()
        
        if Pearson:
            print(">>>SC3 cal Pearson dists")
        
        if Spearman:
            print(">>>SC3 cal Spearman dists")


    def sc3_calc_transformation(self):
        print(">>>SC3 calc transformation")

    def sc3_kmeans(self):
        print(">>>SC3 kmeans cluster")

    def sc3_calc_consensus(self):
        print(">>>SC3 calc consensus matrix")

    def sc3_run_svm(self):
        print(">>>SC3 run SVM (mixed model)") 

    def cal_metric_adjustedRandIndex(self):
        print(">>>SC3 cal metric adjustedRandIndex")


if __name__ == "__main__":
    # os.getcwd()
    print(">>>SC3 algorithm for single cell RNA seq data")

    data_root = "SC3/Data/Yan/yan_export_from_R.csv"
    anno_root = "SC3/Data/Yan/cell_types_export_from_R.txt"
    sc3 = SC3(data_root, anno_root=anno_root)
    
    sc3.sc3_preprocess()
    est_k = sc3.sc3_estimate_k()
    