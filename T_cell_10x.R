library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)

# Load the t_cell dataset
t_cell.data <- Read10X(data.dir = "C:/Users/alexj/Documents/UCD/Research Project/Data Sets/filtered_gene_bc_matrices_t_cell/GRCh38")
# Initialize the Seurat object with the raw (non-normalized data).
t_cell <- CreateSeuratObject(counts = t_cell.data, project = "t_cell3k", min.cells = 3, min.features = 200)
#?CreateSeuratObject
t_cell

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
t_cell[["percent.mt"]] <- PercentageFeatureSet(t_cell, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(t_cell, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(t_cell, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(t_cell, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") # Linear Relationship observed between cell counts and feature counts 
plot1 + plot2

t_cell <- subset(t_cell, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalise the data
t_cell <- NormalizeData(t_cell, normalization.method = "LogNormalize", scale.factor = 10000)

# Highly variable features 
t_cell <- FindVariableFeatures(t_cell, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(t_cell), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(t_cell)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# ?LabelPoints
plot1
plot2

# Scaling the Data 
all.genes <- rownames(t_cell)
t_cell <- ScaleData(t_cell, features = all.genes)

# Linear dimension reduction 
t_cell <- RunPCA(t_cell, features = VariableFeatures(object = t_cell))

# Examine and visualize PCA results a few different ways
print(t_cell[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(t_cell, dims = 1:2, reduction = "pca")

DimPlot(t_cell, reduction = "pca")

DimHeatmap(t_cell, dims = 1:15, cells = 500, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
t_cell <- JackStraw(t_cell, num.replicate = 100)
t_cell <- ScoreJackStraw(t_cell, dims = 1:20)
JackStrawPlot(t_cell, dims = 1:20)

ElbowPlot(t_cell,ndims = 30)
#?ElbowPlot

# Cluster the Cells

t_cell <- FindNeighbors(t_cell, dims = 1:20)
t_cell <- FindClusters(t_cell, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(t_cell), 5)

#install.packages("reticulate")
#reticulate::py_install(packages = 'umap-learn')

t_cell <- RunUMAP(t_cell, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(t_cell, reduction = "umap")

#saveRDS(t_cell, file = )

# Find differently expresses features (cluster Biomarkers)

# find all markers of cluster 2
cluster2.markers <- FindMarkers(t_cell, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(t_cell, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
t_cell.markers <- FindAllMarkers(t_cell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
t_cell.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

cluster0.markers <- FindMarkers(t_cell, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(t_cell, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(t_cell, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(t_cell, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

t_cell.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(t_cell, features = top10$gene) + NoLegend()

