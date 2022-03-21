import sys
import os

import time
import math

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn import preprocessing, decomposition, svm
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
    """
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5410170/
    """
    # init func. read data, anno and determine the n_dims for subsequent transformation components
    @runtime_statistics
    def __init__(self, data_root=None, anno_root=None, d_region_min=0.04, d_region_max = 0.07, MAX_DIM=15, 
                pct_dropout_min=0.1, pct_dropout_max=0.9, kmeans_nstart=1000, kmeans_iter_max=1e+09, 
                times=10, k_range=None, num_SVM_train=0):
        """
            @Parameters  \n
            data_root: full data root for gene-cell scRNAseq data file [default: None]\n
            anno_root: full anno root for corresponding `cell_types` file of gene-cell scRNAseq data [default: None]\n   
            d_region_min: for final feats dims range. min ratio of cell samples [default: 0.04]\n
            d_region_max: for final feats dims range. max ratio of cell samples [default: 0.07]\n
            MAX_DIM: for final feats dims MAX limitation. [default: 15]\n
            pct_dropout_min: the min percent for inital gene filter [default: 0.1]\n
            pct_dropout_max: the max percent for inital gene filter [default: 0.9]\n
            kmeans_nstart: k-means cluster `RANDOMSTATE`param. [default: 50]\n
            kmeans_iter_max: k-means cluster max iteration limitation [default: 1e+09]\n
            times: k-means `n_init` param. Number of time k-means will be run with different centroid seeds [default: 10]\n
            k_range: the range of the number of `k` centroids for kmeans clustering, 
                                            if not set, call `sc3_estimate_k()` for estimation k [default: None]\n
            num_SVM_train: number of samples for SVM training on SVM-mixed workflow [default: 0]\n
        """
        # print(">>>SC3 init")
        assert data_root is not None, "please pass data root!"
        self.adata = sc.read_csv(data_root).transpose() # csv文件中行为基因列为cell sc读进来行列倒置了因此需要转置
        self.adata.uns = dict()
        # self.adata.X # matrix val
        # self.adata.obs # cell name, panda dataframe
        # self.adata.var # gene name, panda dataframe

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
        
        self.d_region_min = d_region_min
        self.d_region_max = d_region_max
        self.MAX_DIM = MAX_DIM

        self.pct_dropout_min = pct_dropout_min
        self.pct_dropout_max = pct_dropout_max

        self.kmeans_nstart = kmeans_nstart
        self.kmeans_iter_max = kmeans_iter_max
        self.times = times

        self.k_range = k_range
        
        assert num_SVM_train < self.adata.n_obs and num_SVM_train >= 0, "num_SVM_train should be in [0, adata.n_obs]!"
        self.num_samples4cluster = 5000 if self.adata.n_obs > 5000 else int(num_SVM_train)
        if self.num_samples4cluster == 0:
            self.num_samples4cluster = self.adata.n_obs

        self.sc3_prepare_cluster_flag()
        
        self.sc3_determine_n_dims()

    
    # one go run 
    @runtime_statistics
    def sc3_onego_run(self):
        # print(">>>SC3 one-go run workflow")
        self.sc3_run_cluster_workflow()
        self.sc3_run_svm()
        self.cal_metric_ARI_global()


    # run the entire clustering workflow, regardless of svm training.
    @runtime_statistics
    def sc3_run_cluster_workflow(self):
        # print(">>>SC3 run clustering workflow")
        self.sc3_preprocess()
        self.sc3_calc_dists()
        self.sc3_calc_transformation()
        self.sc3_kmeans()
        self.sc3_calc_consensus()
        self.cal_metric_ARI_only_cluster()
    

    # define dimensions according to number of cells
    @runtime_statistics
    def sc3_determine_n_dims(self):
        num_cell = self.num_samples4cluster

        min_dim = np.floor(self.d_region_min * num_cell).astype(np.uint8)
        max_dim = np.ceil(self.d_region_max * num_cell).astype(np.uint8)
        n_dim = range(min_dim, max_dim + 1) # calc range of dims
        
        # for large datasets restrict the region of dimensions to 15
        if len(n_dim) > self.MAX_DIM:
            n_dim = np.random.choice(n_dim, self.MAX_DIM, replace=False)# without replacement
        
        self.adata.uns["n_dims"] = n_dim


    # self.adata has been gene filter and log2(X+1) trans. has been aligned to R code.
    @runtime_statistics
    def sc3_preprocess(self):
        # print(">>>SC3 preprocessing")
        num_cell = self.adata.shape[0]

        sc.pp.filter_genes(self.adata, min_cells=self.pct_dropout_min*num_cell+1) # align to R code, not equal 
        sc.pp.filter_genes(self.adata, max_cells=self.pct_dropout_max*num_cell-1) # gene filter between (min, max)
        sc.pp.log1p(self.adata, base=2) # log norm
        # print(">>>Matrix shape:", self.adata.shape)


    # Estimate the optimal k for k-means clustering: has been aligned to R code
    @runtime_statistics
    def sc3_estimate_k(self):
        # print(">>>SC3 estimate k val")
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
        # print(">>>SC3 cal dists")
        data_matrix = self.adata.X[self.adata.obs["samples_cluster_flag"], :]

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
        # print(">>>SC3 calc transformation")
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
    def sc3_kmeans(self):
        # print(">>>SC3 kmeans cluster")
        if self.k_range is None:
            self.sc3_estimate_k()
            self.k_range = [self.adata.uns["est_k"]]
        
        elif isinstance(self.k_range, (list, tuple)):
            pass
        else:
            try:
                self.k_range = list(self.k_range)
            except:
                print("Invalid k_range in `sc3_kmeans` func, should be `None, list or tuple or iterator`")

        trans_matrix = self.adata.uns["trans_matrix"]
        
        all_ks_sim_matrix = {}
        for k_val in self.k_range:  # 多个不同聚类k
            kmeans_cls = KMeans(n_clusters=k_val, max_iter=int(self.kmeans_iter_max), n_init=self.times, random_state=self.kmeans_nstart)
            
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


    # calculate the consensus matrix on similarity matrix by hcluster method
    @runtime_statistics
    def sc3_calc_consensus(self):
        # print(">>>SC3 calc consensus matrix")
        all_ks_sim_matrix = self.adata.uns["all_ks_sim_matrix"]

        hcluster_res = {}
        for key_, matrix in all_ks_sim_matrix.items():
            hcluster = AgglomerativeClustering(int(key_), linkage="complete")
            hcluster.fit(matrix)
            hcluster_res[key_] = hcluster.labels_
        self.adata.uns["hcluster_res"] = hcluster_res


    # prepare cluster data flag, svm-mixed or not supported
    @runtime_statistics
    def sc3_prepare_cluster_flag(self):
        """
        @note: if enable svm-mixed model train, then `samples_cluster_flag` is flag on k-means cluster(True) or svm pred(False).
                if not,  `samples_cluster_flag` is all-true for k-means cluster.
        """
        # print(">>>SC3 prepare cluster flag for workflow, support SVM-mixed flag set")

        samples_cluster_flag = np.zeros(self.adata.n_obs, dtype=np.bool)
        index = np.random.choice(range(self.adata.n_obs), self.num_samples4cluster, replace=False)
        samples_cluster_flag[index] = True
    
        self.adata.obs["samples_cluster_flag"] = samples_cluster_flag


    # run SVM model on clustered `self.num_samples4cluster` samples for training and predict the remaining.
    @runtime_statistics
    def sc3_run_svm(self):
        if self.num_samples4cluster == self.adata.n_obs: # 仅使用svm-mixed model时才进行run svm
            return
        
        # print(">>>SC3 run svm model based on partial data and pseudo clustered label")
        svm_trainer = svm.LinearSVC()

        part_clustered_res = self.adata.uns["hcluster_res"]
        global_res = {}
        for key_k_val, pseudo_label in part_clustered_res.items():
            tmp_sample_idx = np.array(range(self.adata.n_obs)) # 缓存所有样本的索引
            tmp_global_out = np.zeros(self.adata.n_obs, dtype=np.uint8) # 缓存一个所有样本输出结果数组
            tmp_global_out[tmp_sample_idx[self.adata.obs["samples_cluster_flag"]]] = pseudo_label #先把聚类输出按位置填进去

            train_data = self.adata.X[self.adata.obs["samples_cluster_flag"]]
            test_data  = self.adata.X[~self.adata.obs["samples_cluster_flag"]]
            
            svm_trainer.fit(train_data, pseudo_label)
            pred_cls = svm_trainer.predict(test_data)
            tmp_global_out[tmp_sample_idx[~self.adata.obs["samples_cluster_flag"]]] = pred_cls #再把svm输出按位置填进去
            
            global_res[key_k_val] = tmp_global_out
        
        self.adata.uns["global_res"] = global_res


    # calculate ARI metric globally, 不管是否使用svm-mixed model
    @runtime_statistics
    def cal_metric_ARI_global(self):
        # print(">>>SC3 cal metric adjustedRandIndex globally")
        if self.num_samples4cluster == self.adata.n_obs: # 仅使用svm-mixed model时才计算global ARI指标
            return
        
        labels = self.adata.obs["category"]
        
        global_res = self.adata.uns["global_res"]
        ARI = {}
        for key_, pred in global_res.items():
            ARI[key_] = metrics.adjusted_rand_score(labels.tolist(), pred)
        print(ARI)      


    # calculate ARI metric only on cluster res
    @runtime_statistics
    def cal_metric_ARI_only_cluster(self):
        # print(">>>SC3 cal metric adjustedRandIndex only on cluster results")
    
        labels = self.adata.obs["category"][self.adata.obs["samples_cluster_flag"]]
        
        ks_res = self.adata.uns["hcluster_res"]
        ARI = {}
        for key_, pred in ks_res.items():
            ARI[key_] = metrics.adjusted_rand_score(labels.tolist(), pred)
        print(ARI)        


@runtime_statistics
def main(data_root, anno_root, Krange=None, num4SVM=0):
    # os.getcwd()
    print(">>>SC3 algorithm for single cell RNA seq data")

    sc3 = SC3(data_root, anno_root, k_range=Krange, num_SVM_train=num4SVM)
    sc3.sc3_onego_run()


if __name__ == "__main__":
    np.random.seed(2022)
    
    root = "SC3/Data/"    
    data_root_list = ["Yan/yan_export_from_R.csv", "biase/biase_export_from_R.csv"]
    anno_root_list = ["Yan/cell_types_export_from_R.txt", "biase/cell_types_export_from_R.txt"]

    data_root = [os.path.join(root, filename) for filename in data_root_list]
    anno_root = [os.path.join(root, filename) for filename in anno_root_list]

    idx = 1
    main(data_root[idx], anno_root[idx], Krange=None, num4SVM=40)