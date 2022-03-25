# library()
# install.packages("languageserver")
# install.packages("mclust")
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# install.packages("Rcpp")
# install.packages("RcppArmadillo")
# BiocManager::install("scater")

## ----knitr-options, echo=FALSE, message=FALSE, warning=FALSE------------------
library(knitr)
set.seed(2022)
opts_chunk$set(fig.align = "center", fig.width = 6, fig.height = 5, dev = "png")

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(SingleCellExperiment)
library(SC3)
library(scater)

ONEGO <- FALSE
Alldata <- c("Yan", "Biase", "Goolam", "Deng", "Baron")
DATA_NAME <- Alldata[4]
if (DATA_NAME == "Yan") {
    head(ann)
    yan[1:3, 1:3] # vis yan data
    # write.table(ann, "SC3/Data/Yan/cell_types_export_from_R.txt")
    # write.csv(yan, "SC3/Data/Yan/yan_export_from_R.csv")
    # create a SingleCellExperiment object
    sce <- SingleCellExperiment(
        assays = list(
            counts = as.matrix(yan),
            logcounts = log2(as.matrix(yan) + 1)
        ),
        colData = ann
    )
    # define feature names in feature_symbol column
    rowData(sce)$feature_symbol <- rownames(sce)
    # remove features with duplicated names
    sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
} else if (DATA_NAME == "Biase"){
    sce <- readRDS("SC3/Data/biase/biase.rds")
    rownames(sce)
    rowData(sce)
    colnames(sce)
    colData(sce)
    colData(sce)$cell_type1
    colData(sce)$cell_type2
    assayNames(sce)
    counts(sce)[1, ]
    # normcounts(sce)[1, ]
    # logcounts(sce)[1, ]
} else if (DATA_NAME == "Goolam"){
    sce <- readRDS("SC3/Data/Goolam/goolam.rds")
    rownames(sce)
    rowData(sce)
    colnames(sce)
    colData(sce)
    colData(sce)$cell_type1
    colData(sce)$cell_type2
    assayNames(sce)
    # counts(sce)[1, ]
    # logcounts(sce)[1, ]
} else if(DATA_NAME == "Deng"){
    sce <- readRDS("SC3/Data/Deng/deng-rpkms.rds")
    rownames(sce)
    rowData(sce)
    colnames(sce)
    colData(sce)
    colData(sce)$cell_type1
    colData(sce)$cell_type2
    assayNames(sce)
    # logcounts(sce)[1, ]
} else if(DATA_NAME == "Baron"){
    sce <- readRDS("SC3/Data/Baron/baron-human.rds")
    rownames(sce)
    rowData(sce)
    colnames(sce)
    colData(sce)
    colData(sce)$cell_type1
    colData(sce)$cell_type2
    assayNames(sce)
}


if (ONEGO) {
    ## -----------------------------------------------------------------------------
    sce <- runPCA(sce)
    plotPCA(sce, colour_by = "cell_type2")

    ## -----------------------------------------------------------------------------
    sce <- sc3(sce, n_cores = 1, biology = FALSE)

    ## ---- eval=FALSE--------------------------------------------------------------
    sc3_interactive(sce)

    ## ----eval=FALSE---------------------------------------------------------------
    #  sc3_export_results_xls(sce)

    ## -----------------------------------------------------------------------------
    col_data <- colData(sce)
    head(col_data[, grep("sc3_", colnames(col_data))])

    cluter_name <- paste("sc3", metadata(sce)$sc3$k_estimation, "clusters", sep="_")
    ## -----------------------------------------------------------------------------
    sce <- runPCA(sce)
    plotPCA(
        sce,
        colour_by = cluter_name
    )
    # size_by = "sc3_3_log2_outlier_score"
    # )

    ## -----------------------------------------------------------------------------
    row_data <- rowData(sce)
    head(row_data[, grep("sc3_", colnames(row_data))])

    ## ---- fig.height=6------------------------------------------------------------
    sc3_plot_consensus(sce, k = 3)

    ## ---- fig.height=6, fig.width=8-----------------------------------------------
    sc3_plot_consensus(
        sce,
        k = 3,
        show_pdata = c(
            "cell_type1",
            "log10_total_features",
            cluter_name,
            "sc3_3_log2_outlier_score"
        )
    )

    ## -----------------------------------------------------------------------------
    sc3_plot_silhouette(sce, k = 3)

    ## ---- fig.height=6------------------------------------------------------------
    sc3_plot_expression(sce, k = 3)

    ## ---- fig.height=6, fig.width=8-----------------------------------------------
    # sc3_plot_expression(
    #     sce, k = 3,
    #     show_pdata = c(
    #         "cell_type1",
    #         "log10_total_features",
    #         cluter_name,
    #         "sc3_3_log2_outlier_score"
    #     )
    # )

    ## ---- fig.height=3------------------------------------------------------------
    sc3_plot_cluster_stability(sce, k = 3)

    ## ---- fig.height=9------------------------------------------------------------
    # sc3_plot_de_genes(sce, k = 3)

    ## ---- fig.height=9, fig.width=8-----------------------------------------------
    # sc3_plot_de_genes(
    #     sce, k = 3,
    #     show_pdata = c(
    #         "cell_type1",
    #         "log10_total_features",
    #         cluter_name,
    #         "sc3_3_log2_outlier_score"
    #     )
    # )

    ## ---- fig.height=6------------------------------------------------------------
    # sc3_plot_markers(sce, k = 3)

    ## ---- fig.height=6, fig.width=8-----------------------------------------------
    # sc3_plot_markers(
    #     sce, k = 3,
    #     show_pdata = c(
    #         "cell_type1",
    #         "log10_total_features",
    #         cluter_name,
    #         "sc3_3_log2_outlier_score"
    #     )
    # )
} else {
    ## -----------------------------------------------------------------------------
    sce <- sc3_prepare(sce, n_cores = 1, rand_seed = 2022)
    str(metadata(sce)$sc3)

    ## -----------------------------------------------------------------------------
    sce <- sc3_estimate_k(sce)
    str(metadata(sce)$sc3)
    est_k <- metadata(sce)$sc3$k_estimation
    ks_val <- c(metadata(sce)$sc3$k_estimation)

    ## -----------------------------------------------------------------------------
    sce <- sc3_calc_dists(sce)
    names(metadata(sce)$sc3$distances)

    ## -----------------------------------------------------------------------------
    sce <- sc3_calc_transfs(sce)
    names(metadata(sce)$sc3$transformations)
    
    if (FALSE){
        ## Rcpp call C norm_laplacian code
        library(Rcpp)
        library(RcppArmadillo)
        sourceCpp(file='SC3/R/norm_laplacian.cpp')

        dists <- metadata(sce)$sc3$distances
        eu <- get("euclidean", dists)
        L <- norm_laplacian(eu)
        L[1, 1:10]
        L[2, 1:10]
        l <- eigen(L)
        tmp <- l$vectors[, order(l$values)]

        n_dim <- metadata(sce)$sc3$n_dim

        trans <- tmp[, 1:max(n_dim)]
    }

    ## -----------------------------------------------------------------------------
    metadata(sce)$sc3$distances

    ## -----------------------------------------------------------------------------
    sce <- sc3_kmeans(sce, ks = ks_val)
    # write.csv(metadata(sce)$sc3$kmeans, "tmp.csv")
    names(metadata(sce)$sc3$kmeans)

    ## -----------------------------------------------------------------------------
    col_data <- colData(sce)
    head(col_data[, grep("sc3_", colnames(col_data))])

    ## -----------------------------------------------------------------------------
    sce <- sc3_calc_consens(sce)
    names(metadata(sce)$sc3$consensus)

    names(metadata(sce)$sc3$consensus$est_k)

    ## -----------------------------------------------------------------------------
    metadata(sce)$sc3$kmeans

    ## -----------------------------------------------------------------------------
    col_data <- colData(sce)
    head(col_data[, grep("sc3_", colnames(col_data))])

    ## -----------------------------------------------------------------------------
    # sce <- sc3_calc_biology(sce, ks = 2:4)

    ## -----------------------------------------------------------------------------
    col_data <- colData(sce)
    head(col_data[, grep("sc3_", colnames(col_data))])

    ## -----------------------------------------------------------------------------
    row_data <- rowData(sce)
    head(row_data[, grep("sc3_", colnames(row_data))])

    ## -----------------------------------------------------------------------------
    cluter_name <- paste("sc3", metadata(sce)$sc3$k_estimation, "clusters", sep="_")
    no_svm_labels <- colData(sce)[cluter_name][,1]

    ## -----------------------------------------------------------------------------
    sce <- sc3(sce, ks=ks_val,  n_cores = 1, biology = FALSE, svm_num_cells = 50)

    ## -----------------------------------------------------------------------------
    col_data <- colData(sce)
    head(col_data[, grep("sc3_", colnames(col_data))])

    ## ---- message=FALSE, warning=FALSE--------------------------------------------
    sce <- sc3_run_svm(sce, ks = ks_val)
    col_data <- colData(sce)
    head(col_data[, grep("sc3_", colnames(col_data))])

    ## -----------------------------------------------------------------------------
    metadata(sce)$sc3$svm_train_inds <- NULL
    # sce <- sc3_calc_biology(sce, ks = 2:4)
    col_data <- colData(sce)
    head(col_data[, grep("sc3_", colnames(col_data))])

    ## -----------------------------------------------------------------------------
    svm_labels <- colData(sce)[cluter_name][, 1]

    ## -----------------------------------------------------------------------------
    if (require("mclust")) {
        print(paste("no svm vs. svm:", adjustedRandIndex(no_svm_labels, svm_labels)))
        print(paste("svm vs. label:", adjustedRandIndex(svm_labels, colData(sce)$cell_type1)))
        print(paste("no svm vs. label1:", adjustedRandIndex(no_svm_labels, colData(sce)$cell_type1)))
        print(paste("no svm vs. label2:", adjustedRandIndex(no_svm_labels, colData(sce)$cell_type2)))
    }
}

