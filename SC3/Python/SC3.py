import sys
import os

import time
import math

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn import preprocessing, decomposition
from sklearn import metrics
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.cluster import KMeans, AgglomerativeClustering, Birch

from utils import *


# @summary: cal the run time of the function
def runtime_statistics(func):
    def wrapper(*args,**kw):
        start_time = time.time()
        res = func(*args,**kw)
        end_time = time.time()
        print('RUN TIME [{}]: {:.4f} s'.format(func.__name__, (end_time - start_time)))
        return res
    return wrapper


# SC3 algorithm for single cell RNA seq data
class SC3(object):
    # init func. read data, anno and determine the n_dims for subsequent transformation components
    @runtime_statistics
    def __init__(self, data_root=None, anno_root=None):
        print(">>>SC3 init")
        assert data_root is not None, "please pass data root!"
        self.adata = sc.read_csv(data_root).transpose() # csv文件中行为基因列为cell sc读进来行列倒置了因此需要转置
        self.adata.uns = dict()
        # self.adata.X # matrix val
        # self.adata.obs # cell name, panda dataframe
        # self.adata.var # gene name, panda dataframe
        self.sc3_determine_n_dims()

        if anno_root is not None:
            # with open(anno_root, "r") as f:
            #     anno = f.readlines()
            # anno = dict([(each.split("\"")[1], each.split("\"")[3]) for each in anno[1:]]) # {gene_name: cell_type} dict
            anno = pd.read_csv(anno_root, sep=" ")

            # self.adata.uns["anno"] = anno
            if "biase" in anno_root:
                del anno["cell_type2"]

            self.adata.obs["cell_type"] = anno # set to obs
            self.adata.uns["type_names"] = anno["cell_type1"].unique().tolist() # save type names
            self.adata.obs["category"] = anno["cell_type1"].apply(lambda x: self.adata.uns["type_names"].index(x)) # categorial

        else:
            # anno = []
            anno = pd.DataFrame()       
    
    @runtime_statistics
    def sc3_one_go_run(self, k_range=None):
        print(">>>SC3 one-go running")
        self.sc3_preprocess()
        self.sc3_estimate_k()
        self.sc3_calc_dists()
        self.sc3_calc_transformation()
        self.sc3_kmeans(k_range=k_range)
        self.sc3_calc_consensus()
        self.cal_metric_adjustedRandIndex()
    
    # define dimensions according to number of cells
    @runtime_statistics
    def sc3_determine_n_dims(self, d_region_min=0.04, d_region_max = 0.07, MAX_DIM=15):
        num_cell = self.adata.shape[0]
        min_dim = np.floor(d_region_min * num_cell).astype(np.uint8)
        max_dim = np.ceil(d_region_max * num_cell).astype(np.uint8)
        n_dim = range(min_dim, max_dim + 1) # calc range of dims
        
        # for large datasets restrict the region of dimensions to 15
        if len(n_dim) > MAX_DIM:
            n_dim = np.random.choice(n_dim, MAX_DIM, replace=False)# without replacement
        
        self.adata.uns["n_dims"] = n_dim

    # self.adata has been gene filter and log2(X+1) trans. has been aligned to R code.
    @runtime_statistics
    def sc3_preprocess(self, pct_dropout_min=0.1, pct_dropout_max=0.9):
        print(">>>SC3 preprocessing")
        num_cell = self.adata.shape[0]

        sc.pp.filter_genes(self.adata, min_cells=pct_dropout_min*num_cell+1) # align to R code, not equal 
        sc.pp.filter_genes(self.adata, max_cells=pct_dropout_max*num_cell-1) # gene filter between (min, max)
        sc.pp.log1p(self.adata, base=2) # log norm
        print(">>>Matrix shape:", self.adata.shape)
        
    # Estimate the optimal k for k-means clustering: has been aligned to R code
    @runtime_statistics
    def sc3_estimate_k(self):
        print(">>>SC3 estimate k val")
        num_cell, num_gene = self.adata.shape
        muTW = pow((math.sqrt(num_gene - 1) + math.sqrt(num_cell)), 2)
        sigmaTW = (math.sqrt(num_gene - 1) + math.sqrt(num_cell)) * pow((1/math.sqrt(num_gene - 1) + 1/math.sqrt(num_cell)), 1/3)
        bd = 3.273 * sigmaTW + muTW  # 3.2730 is the p=0.001 percentile point for the Tracy-Widom distribution

        x = preprocessing.scale(self.adata.X, axis=1) # warning TB fixed 
        sigmaHatNaive = np.matmul(x, x.transpose())
        
        # compute eigenvalues and return the amount which falls above the bound
        # evals, _ = np.linalg.eig(sigmaHatNaive)
        _, evals, _ = np.linalg.svd(sigmaHatNaive)
        # k = 0
        # for i in range(len(evals)):
        #     if evals[i] > bd:
        #         k += 1
        k = (evals > bd).sum()
        self.adata.uns["est_k"] = k

    # Calculate a distance matrix: has been aligned to R code
    #  Distance between the cells are calculated using `Euclidean, Pearson and Spearman` metrics to construct.
    @runtime_statistics
    def sc3_calc_dists(self, Euclidean=True, Pearson=True, Spearman=True):
        print(">>>SC3 cal dists")
        data_matrix = self.adata.X

        dists = []
        dist_names = []
        if Euclidean:
            print(">>>SC3 cal dists: Euclidean")
            dists.append(pairwise_distances(data_matrix, metric="euclidean"))
            dist_names.append("Euclidean")
        
        if Pearson: # pearson cor 取值范围为 -1 1 因此pearson_distance 取值范围为 0 2
            print(">>>SC3 cal dists: Pearson")
            CENTERED = True # align to R code

            if CENTERED:
                dists.append(pairwise_distances(data_matrix, metric=pearson_distance_centered))
            else:
                dists.append(pairwise_distances(data_matrix, metric=pearson_distance_nocentered))
            dist_names.append("Pearson")

        if Spearman:
            print(">>>SC3 cal dists: Spearman")
            dists.append(pairwise_distances(data_matrix, metric=spearman_distance))
            dist_names.append("Spearman")
        
        if len(dists):
            self.process_data = np.stack(dists, axis=-1) # stack and add new axis 
            self.dist_names = dist_names
        else:
            self.process_data = None
            self.dist_names = None
        assert self.process_data is not None, "Please set dists calculation TYPE!"

    # Calculate transformations of the distance matrices: has been aligned to R code
    @runtime_statistics
    def sc3_calc_transformation(self, pca=True, laplacian=True):
        print(">>>SC3 calc transformation")
        dists_matrix = self.process_data
        n_component = max(self.adata.uns["n_dims"])

        trans_dicts = {}
        if pca:
            print(">>>SC3 calc transformation: PCA")

            pca_cls = decomposition.PCA(n_component)
            # decomposition.KernelPCA(dists_matrix)
            for idx in range(dists_matrix.shape[-1]):
                key_ = "pca_" + self.dist_names[idx]

                pca_cls.fit(preprocessing.scale(dists_matrix[..., idx], axis=0, copy=True)) # has been aligned to R code
                # _ = pca_cls.transform(dists_matrix[..., idx])
                trans_dicts[key_] = pca_cls.components_.transpose()                    

        if laplacian:
            print(">>>SC3 calc transformation: LAPLACIAN")

            for idx in range(dists_matrix.shape[-1]):
                key_ = "laplacian_" + self.dist_names[idx]
                components_ = norm_laplacian(dists_matrix[..., idx], n_component) # has been aligned to R code
                trans_dicts[key_] = components_
        
        self.adata.uns["trans_matrix"] = trans_dicts # trans_dicts are set to anndata.uns for k-means

    # clustering of cells by k-means.
    @runtime_statistics
    def sc3_kmeans(self, k_range=None, times=None, kmeans_nstart = 1000, kmeans_iter_max = 1e+09):
        print(">>>SC3 kmeans cluster")
        if k_range is None:
            k_range = [self.adata.uns["est_k"]]
        elif isinstance(k_range, (list, tuple)):
            pass
        else:
            try:
                k_range = list(k_range)
            except:
                print("Invalid k_range in `sc3_kmeans` func, should be `None, list or tuple or iterator`")
        
        times = 10 if times is None else int(times)

        trans_matrix = self.adata.uns["trans_matrix"]
        
        all_ks_sim_matrix = {}
        for k_val in k_range:  # 多个不同聚类k
            kmeans_cls = KMeans(n_clusters=k_val, max_iter=int(kmeans_iter_max), n_init=times)
            
            k_cluster_matrix = {}
            for key_, val_matrix in trans_matrix.items(): # 单个k下的dist+transformation组合遍历
                for d_dim in self.adata.uns["n_dims"]: #单个k下的某一dist+transformation组合的前d_dim维
                    kmeans_cls.fit(val_matrix[:, :d_dim])

                    matrix_key = "_".join(["d", str(d_dim), key_])
                    k_cluster_matrix[matrix_key] = kmeans_cls.labels_
            
            assert len(k_cluster_matrix.keys()) == len(self.adata.uns["n_dims"])*len(trans_matrix.keys()) 
            k_sim_matrix = convert_similarity_matrix(k_cluster_matrix, n_cells=val_matrix.shape[0], n_clusters=k_val)
            all_ks_sim_matrix[str(k_val)] = k_sim_matrix

        self.adata.uns["all_ks_sim_matrix"] = all_ks_sim_matrix

    @runtime_statistics
    def sc3_calc_consensus(self):
        print(">>>SC3 calc consensus matrix")
        all_ks_sim_matrix = self.adata.uns["all_ks_sim_matrix"]

        hcluster_res = {}
        for key_, matrix in all_ks_sim_matrix.items():
            hcluster = AgglomerativeClustering(int(key_), linkage="complete")
            hcluster.fit(matrix)
            hcluster_res[key_] = hcluster.labels_
        self.adata.uns["hcluster_res"] = hcluster_res

    @runtime_statistics
    def sc3_run_svm(self):
        print(">>>SC3 run SVM (mixed model)")
        pass

    @runtime_statistics
    def cal_metric_adjustedRandIndex(self):
        print(">>>SC3 cal metric adjustedRandIndex")

        labels = self.adata.obs["category"]
        
        ks_res = self.adata.uns["hcluster_res"]
        ARI = {}
        for key_, pred in ks_res.items():
            ARI[key_] = metrics.adjusted_rand_score(labels.tolist(), pred)
        print(ARI)        


@runtime_statistics
def main(data_root, anno_root=None, ONEGO=False, K=None):
    # os.getcwd()
    print(">>>SC3 algorithm for single cell RNA seq data")

    sc3 = SC3(data_root, anno_root=anno_root)
    
    if ONEGO:
        sc3.sc3_one_go_run(k_range=K)
    else:
        sc3.sc3_preprocess()
        sc3.sc3_estimate_k()
        sc3.sc3_calc_dists()
        sc3.sc3_calc_transformation()
        sc3.sc3_kmeans(k_range=K)
        sc3.sc3_calc_consensus()
        sc3.cal_metric_adjustedRandIndex()


if __name__ == "__main__":
    np.random.seed(2022)
    
    root = "SC3/Data/"    
    data_root_list = ["Yan/yan_export_from_R.csv", "biase/biase_export_from_R.csv"]
    anno_root_list = ["Yan/cell_types_export_from_R.txt", "biase/cell_types_export_from_R.txt"]

    data_root = [os.path.join(root, filename) for filename in data_root_list]
    anno_root = [os.path.join(root, filename) for filename in anno_root_list]

    idx = 1
    main(data_root[idx], anno_root[idx], ONEGO=True, K=None)
