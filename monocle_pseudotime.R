library(monocle)
library(cellrangerRkit)

# load expression matrix through the load_cellranger_matrix function from cellrangerRkit
cellranger_pipestance_path <- "path/to/analysis/folder"
gbm <- load_cellranger_matrix(cellranger_pipestance_path)

# create a cell dataset with the count data as expression values. 
# Set the lower detection limit to 0.5 -> include all counts.
# Use a negative binomial distribution with fixed variance as expression family, as recommended for UMI data in the monocle documentation.
gbm_cds <- newCellDataSet(exprs(gbm),
                          phenoData = new("AnnotatedDataFrame", data = pData(gbm)),
                          featureData = new("AnnotatedDataFrame", data = fData(gbm)),
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())

# count genes with an average expression of at least 0.1
gbm_cds <- detectGenes(gbm_cds, min_expr = 0.1)

# estimate size factors and dispersions for clustering and dimension reduction
gbm_cds <- estimateSizeFactors(gbm_cds)
gbm_cds <- estimateDispersions(gbm_cds)
disp_table <- dispersionTable(gbm_cds)

# filter the genes to be used for clustering to only include those with a mean expression over all cells of at least 0.1
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)

# mark genes to be used for clustering
gbm_cds <- setOrderingFilter(gbm_cds, unsup_clustering_genes$gene_id)

# plot dispersion vs. mean expression
plot_ordering_genes(gbm_cds)

# plot the variance explained by each principal component
plot_pc_variance_explained(gbm_cds, return_all = F)

# calculate dimensionality reduction by t-stochastic neighbor embedding (t-SNE)
# reduce down to two components
# use 6 PCs as input
gbm_cds <- reduceDimension(gbm_cds, max_components = 2, num_dim = 6,
                        reduction_method = 'tSNE', verbose = T)

# cluster cells into 9 clusters and plot them
gbm_cds <- clusterCells(gbm_cds, num_clusters = 10)
plot_cell_clusters(gbm_cds, 1, 2)

# Pseudotime calculation

# select only genes for pseutotime calculation that are expressed in more than 5% of cells
gbm_cds_expressed_genes <-  row.names(subset(fData(gbm_cds),
                                          num_cells_expressed >= 10))

# check genes for differential expression in the clusters
clustering_DEG_genes <- differentialGeneTest(gbm_cds[gbm_cds_expressed_genes,],
                                             fullModelFormulaStr = '~Cluster',
                                             cores = 8)

# for ordering, only keep the 1,000 genes with the highest adjusted p-value
gbm_cds_ordering_genes <-
  row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]

# mark genes for ordering
gbm_cds <- setOrderingFilter(gbm_cds, ordering_genes = gbm_cds_ordering_genes)

# reduce dimensionality through the DDRTree method by Mao et al.
gbm_cds <- reduceDimension(gbm_cds, method = 'DDRTree')

# order cells (reverse = T may not be necessary if ordering is oriented correctly)
gbm_cds <- orderCells(gbm_cds, reverse=T)

# plot the cell trajectory (pseudotime) 
plot_cell_trajectory(gbm_cds, color_by = "Cluster")

# extract pseudotime
pseudotime <- gbm_cds@phenoData@data$Pseudotime
