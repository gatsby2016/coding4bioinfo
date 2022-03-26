import os
import math

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import scanpy as sc
from sklearn import preprocessing, decomposition, svm
from sklearn import metrics
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.cluster import KMeans, AgglomerativeClustering, Birch

from utils import *


# SC3 algorithm for single cell RNA seq data
class SC3(object):
    """
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5410170/
    """
    # init func. read data, anno and determine the n_dims for subsequent transformation components
    @runtime_statistics
    def __init__(self, data_root=None, anno_root=None, d_region_min=0.04, d_region_max = 0.07, MAX_DIM=15, 
                pct_dropout_min=0.1, pct_dropout_max=0.9, kmeans_nstart=2022, kmeans_iter_max=1e+09, 
                times=100, k_range=None, num_SVM_train=0):
        """
            @Parameters  \n
            data_root: full data root for gene-cell scRNAseq data file [default: None]\n
            anno_root: full anno root for corresponding `cell_types` file of gene-cell scRNAseq data [default: None]\n   
            d_region_min: for final feats dims range. min ratio of cell samples [default: 0.04]\n
            d_region_max: for final feats dims range. max ratio of cell samples [default: 0.07]\n
            MAX_DIM: for final feats dims MAX limitation. [default: 15]\n
            pct_dropout_min: the min percent for inital gene filter [default: 0.1]\n
            pct_dropout_max: the max percent for inital gene filter [default: 0.9]\n
            kmeans_nstart: k-means cluster `RANDOMSTATE`param. [default: 2022]\n
            kmeans_iter_max: k-means cluster max iteration limitation [default: 1e+09]\n
            times: k-means `n_init` param. Number of time k-means will be run with different centroid seeds [default: 10]\n
            k_range: the range of the number of `k` centroids for kmeans clustering, 
                                            if not set, call `sc3_estimate_k()` for estimation k [default: None]\n
            num_SVM_train: number of samples for SVM training on SVM-mixed workflow [default: 0]\n
        """
        # print(">>>SC3 init")
        assert data_root is not None, "please pass data root!"
        self.data_root = data_root

        self.adata = sc.read_csv(data_root).transpose() # csv文件中行为基因列为cell sc读进来行列倒置了因此需要转置
        print("Loading data shape: {}".format(self.adata.X.shape))

        duplicated_symbol = self.adata.var_names.duplicated() #针对性（deng数据）增加特征维度去重
        tmp_var = self.adata.var
        self.adata = sc.AnnData(self.adata.X[:, ~duplicated_symbol], 
                                    obs=self.adata.obs, var=tmp_var[~duplicated_symbol])
        print("After duplicating, data shape: {}".format(self.adata.X.shape))

        self.gene_filter_flag = self.gene_filter(self.adata.X, pct_dropout_min, pct_dropout_max) #基于原始输入计算filter_flag

        # NOTE: 经比对 原工作用count处理确定gen_filter_flag, 然后取对应的logcount;其中goolam和Baron中间要加一个CPMnorm
        if "Goolam" in data_root or "Baron" in data_root or "Klein" in data_root: 
            sc.pp.normalize_total(self.adata, target_sum=1e6, exclude_highly_expressed=True) # CPM normalization aligned 

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
            if "Deng" in anno_root:
                # del anno["cell_type2"]
                del anno["cell_type1"]
                anno = anno.rename(columns={"cell_type2":"cell_type1"})
            if "Baron" in anno_root:
                del anno["human"]

            self.adata.obs["cell_type"] = anno # set to obs
            self.adata.uns["type_names"] = anno["cell_type1"].unique().tolist() # save type names
            self.adata.obs["category"] = anno["cell_type1"].apply(lambda x: self.adata.uns["type_names"].index(x)) # categorial

        else:
            # anno = []
            anno = pd.DataFrame()       
        
        self.d_region_min = d_region_min
        self.d_region_max = d_region_max
        self.MAX_DIM = MAX_DIM

        self.kmeans_nstart = kmeans_nstart
        self.kmeans_iter_max = kmeans_iter_max
        self.times = times

        self.k_range = k_range
        
        assert num_SVM_train < self.adata.n_obs and num_SVM_train >= 0, "num_SVM_train should be in [0, adata.n_obs]!"
        self.num_samples4cluster = 5000 if self.adata.n_obs > 5000 else int(num_SVM_train)
        if self.num_samples4cluster == 0:
            self.num_samples4cluster = self.adata.n_obs
        self.adata.uns["num_samples4cluster"] = self.num_samples4cluster

        self.sc3_prepare_cluster_flag()
        
        self.sc3_determine_n_dims()

    
    # one go run 
    @runtime_statistics
    def sc3_onego_run(self, write=False):
        # print(">>>SC3 one-go run workflow")
        self.sc3_run_cluster_workflow()
        self.sc3_run_svm()
        self.cal_metric_ARI_global()
        
        if write:
            self.save_results()


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

        min_dim = np.floor(self.d_region_min * num_cell).astype(np.uint32)
        max_dim = np.ceil(self.d_region_max * num_cell).astype(np.uint32)
        n_dim = list(range(min_dim, max_dim + 1)) # calc range of dims
        
        # for large datasets restrict the region of dimensions to 15
        if len(n_dim) > self.MAX_DIM:
            n_dim = np.random.choice(n_dim, self.MAX_DIM, replace=False)# without replacement
        
        self.adata.uns["n_dims"] = n_dim

    # gene filter or not, return  FLAG 1 is satisfaction and 0 is need to be filtered.
    # @matrix: (numpy.array) the origin input: samples(obs) x gen_dims 
    def gene_filter(self, matrix, pct_dropout_min=0.1, pct_dropout_max=0.9):
        num_cell = matrix.shape[0]
        dropouts = (np.sum(matrix == 0, axis=0) / num_cell)
        gene_filter_flag = (dropouts < pct_dropout_max) & (dropouts > pct_dropout_min)

        return gene_filter_flag # 0 is need to be filtered.


    # self.adata has been gene filter and log2(X+1) trans. has been aligned to R code.
    @runtime_statistics
    def sc3_preprocess(self):
        # print(">>>SC3 preprocessing")
        sc.pp.log1p(self.adata, base=2) # 先进行log norm
        tmp_var = self.adata.var
        self.adata = sc.AnnData(self.adata.X[:, self.gene_filter_flag], var=tmp_var[self.gene_filter_flag], \
                                    obs=self.adata.obs, uns=self.adata.uns) # 然后保留满足条件的gene
        print(">>>After gene_filter, data shape:", self.adata.shape)


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
                # dists.append(pairwise_distances(data_matrix, metric=pearson_distance_centered))
                dists.append(1.0-np.corrcoef(data_matrix))
            else:
                dists.append(pairwise_distances(data_matrix, metric=pearson_distance_nocentered))
            dist_names.append("Pearson")

        if Spearman:
            print(">>>SC3 cal dists: Spearman")
            # dists.append(pairwise_distances(data_matrix, metric=spearman_distance))
            dists.append(matrix_spearman_distance(data_matrix.transpose()))# more efficient
            dist_names.append("Spearman")
        
        if len(dists):
            self.process_data = np.stack(dists, axis=-1) # stack and add new axis 
            self.dist_names = dist_names
        else:
            self.process_data = None
            self.dist_names = None
        assert self.process_data is not None, "Please set dists calculation TYPE!"
        self.adata.uns["dists_matrix"] = self.process_data


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
        elif isinstance(self.k_range, int):
            self.k_range = [self.k_range]
        else:
            try:
                self.k_range = list(self.k_range)
            except:
                print("Invalid k_range in `sc3_kmeans` func, should be `None, list or tuple or iterator`")

        trans_matrix = self.adata.uns["trans_matrix"]
        
        all_ks_sim_matrix = {}
        for k_val in self.k_range:  # 多个不同聚类k
            kmeans_cls = KMeans(n_clusters=k_val, max_iter=int(self.kmeans_iter_max), 
                                                    n_init=self.times, random_state=self.kmeans_nstart)
            
            k_cluster_matrix = {}
            for d_dim in self.adata.uns["n_dims"]: #单个k下的某一dist+transformation组合的前d_dim维
                for key_, val_matrix in trans_matrix.items(): # 单个k下的dist+transformation组合遍历
                    kmeans_cls.fit(val_matrix[:, :d_dim])

                    matrix_key = "_".join(["d", str(d_dim), key_])
                    k_cluster_matrix[matrix_key] = kmeans_cls.labels_
            
            assert len(k_cluster_matrix.keys()) == len(self.adata.uns["n_dims"])*len(trans_matrix.keys()) 
            k_sim_matrix = convert_similarity_matrix(k_cluster_matrix, n_cells=val_matrix.shape[0], n_clusters=k_val)
            """debug validation with Rcode kmeans output
            # in R: run to SC3/R/SC3.R#L224; then `write.csv(metadata(sce)$sc3$kmeans, "tmp.csv")`
            # in Python: 
            rres = sc.read_csv("tmp.csv")
            aris = {}
            for idx, keys in enumerate(k_cluster_matrix.keys()):
                aris[keys] = metrics.adjusted_rand_score(k_cluster_matrix[keys], rres.X[:, idx])
                print(aris[keys])

            for idx, keys in enumerate(k_cluster_matrix.keys()): # 将R中结果填到k_cluster_matrix计算后续,一致对齐,证实仅在kmeans有差异
                k_cluster_matrix[keys] = (rres.X[:, idx] -1).astype(np.int32) 
            """
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
            # sns.clustermap(matrix, method="complete") # only if you wan to show hierarchical_clustering map
            # plt.show()
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

        samples_cluster_flag = np.zeros(self.adata.n_obs, dtype=bool)
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
        self.cal_metric_ARI(labels, ks_res=global_res)  


    # calculate ARI metric only on cluster res
    @runtime_statistics
    def cal_metric_ARI_only_cluster(self):
        # print(">>>SC3 cal metric adjustedRandIndex only on cluster results")
    
        labels = self.adata.obs["category"][self.adata.obs["samples_cluster_flag"]]
        
        ks_res = self.adata.uns["hcluster_res"]
        self.cal_metric_ARI(labels, ks_res=ks_res)


    @staticmethod
    def cal_metric_ARI(labels, ks_res):
        ARI = {}
        for key_, pred in ks_res.items():
            ARI[key_] = metrics.adjusted_rand_score(labels.tolist(), pred)
        print(ARI)


    @runtime_statistics
    def save_results(self):
        sc.write(os.path.join("../SC3/Results/", self.data_root.split("/")[-2]), self.adata, ext="h5ad")