install.packages("languageserver")
## ----knitr-options, echo=FALSE, message=FALSE, warning=FALSE------------------
library(knitr)
opts_chunk$set(fig.align = 'center', fig.width = 6, fig.height = 5, dev = 'png')

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(SingleCellExperiment)
library(SC3)
library(scater)

head(ann)
yan[1:3, 1:3]

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
sce <- runPCA(sce)
plotPCA(sce, colour_by = "cell_type1")

## -----------------------------------------------------------------------------
sce <- sc3(sce, ks = 2:4, biology = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  sc3_interactive(sce)

## ----eval=FALSE---------------------------------------------------------------
#  sc3_export_results_xls(sce)

## -----------------------------------------------------------------------------
col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])

## -----------------------------------------------------------------------------
sce <- runPCA(sce)
plotPCA(
    sce, 
    colour_by = "sc3_3_clusters", 
    size_by = "sc3_3_log2_outlier_score"
)

## -----------------------------------------------------------------------------
row_data <- rowData(sce)
head(row_data[ , grep("sc3_", colnames(row_data))])

## ---- fig.height=6------------------------------------------------------------
sc3_plot_consensus(sce, k = 3)

## ---- fig.height=6, fig.width=8-----------------------------------------------
sc3_plot_consensus(
    sce, k = 3, 
    show_pdata = c(
        "cell_type1", 
        "log10_total_features",
        "sc3_3_clusters", 
        "sc3_3_log2_outlier_score"
    )
)

## -----------------------------------------------------------------------------
sc3_plot_silhouette(sce, k = 3)

## ---- fig.height=6------------------------------------------------------------
sc3_plot_expression(sce, k = 3)

## ---- fig.height=6, fig.width=8-----------------------------------------------
sc3_plot_expression(
    sce, k = 3, 
    show_pdata = c(
        "cell_type1", 
        "log10_total_features",
        "sc3_3_clusters", 
        "sc3_3_log2_outlier_score"
    )
)

## ---- fig.height=3------------------------------------------------------------
sc3_plot_cluster_stability(sce, k = 3)

## ---- fig.height=9------------------------------------------------------------
sc3_plot_de_genes(sce, k = 3)

## ---- fig.height=9, fig.width=8-----------------------------------------------
sc3_plot_de_genes(
    sce, k = 3, 
    show_pdata = c(
        "cell_type1", 
        "log10_total_features",
        "sc3_3_clusters", 
        "sc3_3_log2_outlier_score"
    )
)

## ---- fig.height=6------------------------------------------------------------
sc3_plot_markers(sce, k = 3)

## ---- fig.height=6, fig.width=8-----------------------------------------------
sc3_plot_markers(
    sce, k = 3, 
    show_pdata = c(
        "cell_type1", 
        "log10_total_features",
        "sc3_3_clusters", 
        "sc3_3_log2_outlier_score"
    )
)

## -----------------------------------------------------------------------------
sce <- sc3_prepare(sce)
str(metadata(sce)$sc3)

## -----------------------------------------------------------------------------
sce <- sc3_estimate_k(sce)
str(metadata(sce)$sc3)

## -----------------------------------------------------------------------------
sce <- sc3_calc_dists(sce)
names(metadata(sce)$sc3$distances)

## -----------------------------------------------------------------------------
sce <- sc3_calc_transfs(sce)
names(metadata(sce)$sc3$transformations)

## -----------------------------------------------------------------------------
metadata(sce)$sc3$distances

## -----------------------------------------------------------------------------
sce <- sc3_kmeans(sce, ks = 2:4)
names(metadata(sce)$sc3$kmeans)

## -----------------------------------------------------------------------------
col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])

## -----------------------------------------------------------------------------
sce <- sc3_calc_consens(sce)
names(metadata(sce)$sc3$consensus)
names(metadata(sce)$sc3$consensus$`3`)

## -----------------------------------------------------------------------------
metadata(sce)$sc3$kmeans

## -----------------------------------------------------------------------------
col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])

## -----------------------------------------------------------------------------
sce <- sc3_calc_biology(sce, ks = 2:4)

## -----------------------------------------------------------------------------
col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])

## -----------------------------------------------------------------------------
row_data <- rowData(sce)
head(row_data[ , grep("sc3_", colnames(row_data))])

## -----------------------------------------------------------------------------
no_svm_labels <- colData(sce)$sc3_3_clusters

## -----------------------------------------------------------------------------
sce <- sc3(sce, ks = 2:4, biology = TRUE, svm_num_cells = 50)

## -----------------------------------------------------------------------------
col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])

## ---- message=FALSE, warning=FALSE--------------------------------------------
sce <- sc3_run_svm(sce, ks = 2:4)
col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])

## -----------------------------------------------------------------------------
metadata(sce)$sc3$svm_train_inds <- NULL
sce <- sc3_calc_biology(sce, ks = 2:4)
col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])

## -----------------------------------------------------------------------------
svm_labels <- colData(sce)$sc3_3_clusters

## -----------------------------------------------------------------------------
if (require("mclust")) {
  adjustedRandIndex(no_svm_labels, svm_labels)
}

