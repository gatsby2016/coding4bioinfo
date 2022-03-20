### DATA
d <- read.table("SC3/Data/biase/GSE57249_fpkm.txt", header = T)
d <- d[!duplicated(d$ID),]
rownames(d) <- d$ID
d <- d[,2:ncol(d)]
head(d)
### ANNOTATIONS
ann <- read.table("SC3/Data/biase/biase_cell_types.txt", stringsAsFactors = F)

write.table(ann, "SC3/Data/biase/cell_types_export_from_R.txt")
write.csv(d, "SC3/Data/biase/biase_export_from_R.csv")
### SINGLECELLEXPERIMENT
source("SC3/R/utils/create_sce.R")
sceset <- create_sce_from_normcounts(d, ann)
# convert ensembl ids into gene names
# gene symbols will be stored in the feature_symbol column of fData
sceset <- getBMFeatureAnnos(
    sceset, filters="ensembl_gene_id",
    biomart="ensembl", dataset="mmusculus_gene_ensembl")
# remove features with duplicated names
sceset <- sceset[!duplicated(rowData(sceset)$feature_symbol), ]
saveRDS(sceset, "SC3/Data/biase/biase.rds")
