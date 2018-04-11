library(Seurat)
library(ggplot2)
library(viridis)

# load data
data_path <- '/path/to/data'
testis.data <- Read10X(data.dir = data_path)
testis <- CreateSeuratObject(raw.data = testis.data, project = "testis", min.cells = 5)

# assign biological replicate IDs (in this case, by GEM group)
testis@meta.data$mouse[grepl("-1", colnames(testis@raw.data))] <- "mouse 1"
testis@meta.data$mouse[grepl("-2", colnames(testis@raw.data))] <- "mouse 2"

# calculate the proportion of mitochondrial transcripts
mito.genes <- grep(pattern = "^mt-", x = rownames(x = testis@data), value = TRUE)
percent.mito <- Matrix::colSums(testis@raw.data[mito.genes, ])/Matrix::colSums(testis@raw.data)
testis <- AddMetaData(object = testis, metadata = percent.mito, col.name = "percent.mito")

# first QC plot (Violin plot of gene count, UMI count, and proportion of mitochondrial transcripts)
VlnPlot(object = testis, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, group.by = "mouse", point.size.use	= 0.1,
        cols.use=viridis_pal()(2))

# optional built-in plotting
#par(mfrow = c(1, 2))
#GenePlot(object = testis, gene1 = "nUMI", gene2 = "percent.mito")
#GenePlot(object = testis, gene1 = "nUMI", gene2 = "nGene")

# scatter plot of QC factors
p1 <- ggplot(testis@meta.data, aes(x=nUMI, y=percent.mito, color=mouse)) + geom_point(alpha=0.5) +
  scale_color_viridis(discrete=T, guide=F) + xlab("UMI count") + ylab("Mitochondrial transcript levels")
p2 <- ggplot(testis@meta.data, aes(x=nUMI, y=nGene, color=mouse)) + geom_point(alpha=0.5) +
  scale_color_viridis(discrete=T) + xlab("UMI count") + ylab("Gene count")
plot_grid(p1, p2, rel_widths=c(0.8, 1))
