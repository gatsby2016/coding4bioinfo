import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

from SC3 import SC3
from utils import runtime_statistics


@runtime_statistics
def main(data_root, anno_root, Krange=None, num4SVM=0, LOAD=True, WRITEROOT=None):
    # os.getcwd()
    print(">>>SC3 algorithm for single cell RNA seq data")

    DATASET = data_root.split("/")[-2]
    if LOAD and os.path.isfile(os.path.join("SC3/Results", WRITEROOT, DATASET+".h5ad")):
        print("Find h5ad file in {}, will load it and skip running".format(os.path.join("SC3/Results", WRITEROOT)))

        adata = sc.read(os.path.join("../SC3/Results/", WRITEROOT, DATASET), ext="h5ad")
        # print("est_k value: {}".format(adata.uns["est_k"]))

        print("===========Only Cluster results ARI metrics")
        SC3.cal_metric_ARI(adata.obs["category"][adata.obs["samples_cluster_flag"]], adata.uns["hcluster_res"])
        for key_, val in adata.uns["all_ks_sim_matrix"].items():
            fig = sns.clustermap(val, method="complete", cmap="jet")
            fig.savefig("SC3/Results/{}_clustermap_{}.png".format(DATASET, key_), dpi=400)

        if "global_res" in adata.uns.keys():
            print("===========Global results ARI metrics on mixed model with SVM")
            SC3.cal_metric_ARI(adata.obs["category"], adata.uns["global_res"])

    else:
        print("No need to load, or Do NOT find h5ad file in {}, will run directly".format(os.path.join("SC3/Results")))

        sc3 = SC3(data_root, anno_root, k_range=Krange, num_SVM_train=num4SVM)
        sc3.sc3_onego_run(writeroot=WRITEROOT)


if __name__ == "__main__":
    np.random.seed(2022)
    
    root = "SC3/Data/"    
    data_root_list = ["Yan/yan_export_from_R.csv", "biase/biase_export_from_R.csv", 
                        "Goolam/goolam_export_from_R.csv", "Baron/Baron_human_export_from_R.csv",
                        "Deng/deng_export_from_R_rpkms.csv", "Klein/klein_export_from_R.csv"]
    anno_root_list = ["Yan/cell_types_export_from_R.txt", "biase/cell_types_export_from_R.txt", 
                        "Goolam/cell_types_export_from_R.txt", "Baron/human_cell_types_export_from_R.txt",
                        "Deng/cell_types_export_from_R.txt", "Klein/cell_types_export_from_R.txt"]

    data_root = [os.path.join(root, filename) for filename in data_root_list]
    anno_root = [os.path.join(root, filename) for filename in anno_root_list]

    idx = 3
    print("Handling dataset: ", anno_root[idx])

    main(data_root[idx], anno_root[idx], Krange=None, num4SVM=0, LOAD=True, WRITEROOT="rbf")