import time
import numpy as np
import pandas as pd
from scipy.spatial import distance
import scipy.stats as ss


# @summary: cal the run time of the function
def runtime_statistics(func):
    def wrapper(*args,**kw):
        start_time = time.time()
        res = func(*args,**kw)
        end_time = time.time()
        print('RUN TIME [{}]: {:.4f} s'.format(func.__name__, (end_time - start_time)))
        return res
    return wrapper


# //' Graph Laplacian calculation
# //' 
# //' Calculate graph Laplacian of a symmetrix matrix
# }
def norm_laplacian(matrix, n_component=None):
    A = np.exp(-matrix/ matrix.max())
    D_row = pow(sum(A), -0.5).reshape(1, -1)
    A *= D_row
    A *= D_row.transpose()
    res = np.eye(A.shape[-1], A.shape[-1]) - A # has been aligned to R code

    # vectors, evals, _ = np.linalg.svd(res)
    # inds = np.argsort(evals) # 倒序 从小到大排序
    # if n_component is None:
    #     return vectors[:, inds]
    # else:
    #     return vectors[:, inds[:n_component]]
    evals, vectors = np.linalg.eigh(res) # only real symmetric matrix factorized 
    return vectors[:, -n_component:] # fixed! get top-k eignvaules' eignvectors


# 3 pearson distance calculation method: 1st, pearson_distance_centered aligned to R code
def pearson_distance_centered(x, y):
    return distance.correlation(x, y, centered=True)

def pearson_distance_nocentered(x, y):
    return distance.correlation(x, y, centered=False)

def pearson_distance_pandas(x, y):
    X1=pd.Series(x)
    Y1=pd.Series(y)

    pearson_corr = X1.corr(Y1,method="pearson")
    return 1 - pearson_corr


# calc spearman distance: 1 - spearman_corr
def spearman_distance(x, y):
    x1=pd.Series(x)
    y1=pd.Series(y)
    
    spearman_corr = x1.corr(y1,method='spearman')
    return 1 - spearman_corr

# calc spearman distance: 1 - spearman_corr by scipy.stats package
# @matrix: should be (variables x obs)
def matrix_spearman_distance(matrix): 
    coefficient, _ = ss.spearmanr(matrix)
    return 1 - coefficient    
   

# convert kmeans clustering labels to binary similarity matrix
# refer to: https://github.com/hemberg-lab/sc3s/blob/master/sc3s/_cluster.py#L124
def convert_similarity_matrix(dict_object, n_clusters, n_cells):
    # check labels has correct number of clusters and cells
    n_cells_list = [len(x) for x in dict_object.values()]
    assert n_cells_list.count(n_cells) == len(n_cells_list), "number of cells is not consistent with kmeans result"
    
    n_clusters_list = [len(np.unique(x)) for x in dict_object.values()]
    assert n_clusters_list.count(n_clusters) == len(n_clusters_list), "number of clusters not consistent with true num of cluster"

    # initialise empty zero-value binary_matrix array with the correct shape
    binary_sim_matrix = np.zeros((n_cells, n_cells), dtype=np.float32)

    # for each run, we create the block to obtain sim_matrix
    for (key_, cell_labels) in dict_object.items():
        # create a total of 'n_clusters' columns for this iteration
        b = np.zeros((n_cells, n_clusters), dtype=int)

        # annotate the correct columns for each row/cell
        b[range(0, n_cells), cell_labels] = 1

        # this checks that every cell only has one cluster assignment
        assert np.all(np.sum(b, axis=1) == 1), "some cells have multiple cluster assignments"

        # obtain the binary similarity matrix (90, 90)
        sim_matrix = np.matmul(b, b.transpose())

        #sum to the overall binary consensus matrix
        binary_sim_matrix = binary_sim_matrix + sim_matrix
    
    binary_sim_matrix /= len(dict_object.keys()) # average on total runs
    assert binary_sim_matrix.max() <= 1.0, "Max value in similarity matrix should be less than total runs"
    
    return binary_sim_matrix