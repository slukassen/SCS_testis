library(ggplot2)
library(Seurat)

# set working directory to cellranger results/outs/analysis
setwd("path/to/results/outs/analysis")

# load dimension reduction and clustering results from the 
tsne <- read.delim("tsne\\2_components\\projection.csv", sep=",", header=T)
graphclust <- read.delim("clustering\\graphclust\\clusters.csv", sep=",")
k10 <- read.delim("clustering\\kmeans_10_clusters\\clusters.csv", sep=",")
k9 <- read.delim("clustering\\kmeans_9_clusters\\clusters.csv", sep=",")

# group IDs determined by marker gene expression
k9_groups <- c("RS1", "SC2", "CS", "RS2", "ES", "SC1", "Sertoli", "Spg", "Leydig")
groups <- rep("", nrow(k9))
for (i in 1:length(groups)){
  groups[i] <- k9_groups[k9$Cluster[i]]
}

# combine metadata into single data.frame
projection <- as.data.frame(cbind(tsne, graphclust$Cluster, k10$Cluster, k9$Cluster))
colnames(projection) <- c("Barcode", "TSNE.1", "TSNE.2", "graphbased clustering", "K10", "K9")
projection$`graphbased clustering` <- as.factor(projection$`graphbased clustering`)
projection$K10 <- as.factor(projection$K10)
projection$K9 <- as.factor(projection$K9)

# plot the t-SNE results, colored by clustering
p <- ggplot(projection, aes(x=TSNE.1, y=TSNE.2, color=`graphbased clustering`)) + 
  geom_point(alpha=0.5, shape=16) + labs(color="cluster") + 
  theme(axis.text=element_blank(), axis.title=element_blank(),
        axis.ticks = element_blank())
p

q <- ggplot(projection, aes(x=TSNE.1, y=TSNE.2, color=K10)) + 
  geom_point(alpha=0.5, shape=16) + labs(color="cluster") + 
  theme(axis.text=element_blank(), axis.title=element_blank(),
        axis.ticks = element_blank())
q

r <- ggplot(projection, aes(x=TSNE.1, y=TSNE.2, color=groups)) + 
  geom_point(alpha=0.5, shape=16) + labs(color="cluster") + 
  theme(axis.text=element_blank(), axis.title=element_blank(),
        axis.ticks = element_blank())
r

# read expression data from cellranger .mtx output
testis.data <- Read10X(data.dir = "E:\\HG\\Chromium\\Wipa_Repeats_aggr\\Wipa_Repeats_aggr\\outs\\filtered_gene_bc_matrices_mex\\mm10Repeats")

# create Seurat object, filtering genes for those expressed in at least 3 cells
testis <- CreateSeuratObject(raw.data=testis.data, project="Testis", min.cells=3)

# assign the group IDs (from K-means clustering with K=9)
testis@meta.data$group <- groups

# assign the group IDs as cell identity, ordering the factor by termporal succession of cell types and placing somatic cells in the beginning
testis@ident <- factor(testis@meta.data$group, levels=c("Leydig", "Sertoli", "Spg", "SC1", "SC2", "RS1", "RS2", "ES", "CS"))

# normalize the expression data
testis <- NormalizeData(testis)

# scale the expression data
testis <- ScaleData(testis, display.progress = T)

# load marker genes to plot
markers <- read.delim("path/to/marker/genes", header=T)[,1]

# filter gene names to only include those which are annotated in the expression data
markers <- markers[markers %in% row.names(testis@data)]

# dot plot of marker gene expression in the different clusters
sdp <- DotPlot(testis, genes.plot = markers, cols.use = c("grey", "blue"), 
                      x.lab.rot = T, plot.legend = T, dot.scale = 8, do.return = T)
