library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


load("SRA866994_SRS4545966.sparse.RData")

##change rownames(genes)
r <- rownames(sm)
split <- strsplit(r, "_ENS")
new <- unlist(lapply(split, `[[`, 1))
row.names(sm) <- new

head(rownames(sm))

##SeuratObject
scr <- CreateSeuratObject(counts = sm, project = "single-cell-brain", min.cells = 3, min.features = 200)

##Quality control
grep("^mt-",rownames(scr),value = TRUE)
scr[["percent.mt"]] <- PercentageFeatureSet(scr, pattern = "^mt-")

grep("^Rp[ls]", rownames(scr), value = TRUE)
scr[["percent.rbp"]] <- PercentageFeatureSet(scr, pattern = "^Rp[ls]")

VlnPlot(scr, features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rbp"), pt.size = 0, ncol = 4)
plot1 <- FeatureScatter(scr, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scr, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(scr, feature1 = "nCount_RNA", feature2 = "percent.rbp")
plot1
plot2
plot3

dim(scr)

scr <- subset(scr, subset = nFeature_RNA > 250 & nFeature_RNA < 3500 & percent.mt < 5)

dim(scr)

##NORMALIZATION
scr <- NormalizeData(scr, normalization.method = "LogNormalize", scale.factor = 10000)


##ADDING CELL CYCLE PHASES
cc.genes.updated.2019
CellCycleScoring(scr, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE) -> scr


##PCA
scr <- FindVariableFeatures(scr, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(scr), 10)

plot4 <- VariableFeaturePlot(scr)
plot5 <- LabelPoints(plot = plot4, points = top10, repel = TRUE)
plot4
plot5

all.genes <- rownames(scr)
scr <- ScaleData(scr, features = all.genes)
scr <- RunPCA(scr, features = VariableFeatures(object = scr))
print(scr[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(scr, dims = 1:2, reduction = "pca", ncol=2)

##Cell cycle in PC1 AND PC2
DimPlot(scr, dims = 1:2, reduction = "pca")


##Choosing dimensions
ElbowPlot(scr, ndims=40)

pc.touse <- (scr$pca@stdev)^2
pc.touse <- pc.touse/sum(pc.touse)
pc.touse <- cumsum(pc.touse)[1:50]
pc.touse <- min(which(pc.touse>=0.75))
pc.touse


##CLUSTERING

##dimensions till 14
scr_1 <- FindNeighbors(scr, dims = 1:14)
scr_1_0.5 <- FindClusters(scr_1, resolution = 0.5)
scr_1_0.8 <- FindClusters(scr_1, resolution = 0.8)

table(Idents(scr_1_0.5))
head(scr_1_0.5[[]],5)

table(Idents(scr_1_0.8))
head(scr_1_0.8[[]],5)

cell.num <- table(Idents(scr_1_0.8))

ClusterLabels = paste("Cluster", names(cell.num), paste0("(n = ", cell.num, ")"))

ClusterBreaks = names(cell.num)

DimPlot(scr_1_0.8, reduction = "pca")
DimPlot(scr_1_0.8,reduction="pca", dims=c(5,6))

scr_1_0.8 <- RunTSNE(scr_1_0.8, dims=1:14)
DimPlot(scr_1_0.8, reduction = "tsne")

scr_1_0.8 <- RunUMAP(scr_1_0.8, dims = 1:14)
DimPlot(scr_1_0.8, reduction = "umap") + scale_colour_discrete(breaks = ClusterBreaks, labels = ClusterLabels)

##dimensions till 20
scr <- FindNeighbors(scr, dims = 1:20)
scr <- FindClusters(scr, resolution = 0.8)

table(Idents(scr))

cell.num <- table(Idents(scr))

ClusterLabels = paste("Cluster", names(cell.num), paste0("(n = ", cell.num, ")"))

ClusterBreaks = names(cell.num)

DimPlot(scr, reduction = "pca")
DimPlot(scr,reduction="pca", dims=c(5,6))

scr <- RunTSNE(scr, dims=1:20)
DimPlot(scr, reduction = "tsne")

scr <- RunUMAP(scr, dims = 1:20)
DimPlot(scr, reduction = "umap", label = TRUE) + scale_colour_discrete(breaks = ClusterBreaks, labels = ClusterLabels)


##Check biases:
VlnPlot(scr,features="percent.rbp")
VlnPlot(scr,features="nFeature_RNA")
VlnPlot(scr,features="percent.mt")
VlnPlot(scr,features="nCount_RNA")

scr@meta.data %>%
  group_by(seurat_clusters,Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster")


##Identifying marker genes:
scr.markers <- FindAllMarkers(scr, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

scr.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

scr.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(scr, features = top10$gene) + NoLegend()

cluster21.markers <- FindMarkers(scr, ident.1 = 21, min.pct = 0.25, test.use = "wilcox")
cluster21.markers <- cluster21.markers[order(-cluster21.markers$avg_log2FC),]
head(cluster21.markers, n = 10)

cluster1516.markers <- FindMarkers(scr, ident.1 = 15, ident.2 = 16, min.pct = 0.25, test.use = "wilcox")
cluster1516.markers <- cluster1516.markers[order(cluster1516.markers$avg_log2FC),]
head(cluster1516.markers, n = 10)

cluster5AND15.markers <- FindMarkers(scr, ident.1 = c(5,15), min.pct = 0.25, test.use = "wilcox")
cluster5AND15.markers <- cluster5AND15.markers[order(-cluster5AND15.markers$avg_log2FC),]
head(cluster5AND15.markers, n = 10)

FeaturePlot(scr, features = c("C1qa", "Camp", "Ly6c2", "Ms4a1","Cd19","Siglech", "Xcr1", "Cd209a"))

FeaturePlot(scr, features = c("Ccr9","Ccr7", "Ace", "Nkg7", "Klrb1c","Cd3e", "Cd3d", "Adgre1", "Fcgr1"))

FeaturePlot(scr, features = c("Nudt17"))

DotPlot(scr, features = c("C1qa","Camp","Ms4a1","Siglech", "Xcr1", "Cd209a","Cd3e","Ace","Ly6c2","Adgre1", "Fcgr1","Klrb1c"))


##Constructing UMAP with known cells
new.cluster.ids <- c("B", "Macrophages", "Neutrophils", "Macrophages", "Macrophages", "Neutrophils", "NK", "cDC2", "B", "B", "Monocyte", "pDC", "T", "cDC1", "nc monocyte", "Neutrophils", "Monocyte", "T","T","Monocyte", "migDC", "Unknown")
names(new.cluster.ids) <- levels(scr)
scr <- RenameIdents(scr, new.cluster.ids)
DimPlot(scr, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



