# # 安装并导入必要的包
# options("repos" = c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"));options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor");
# if (!"Seurat" %in% installed.packages()) {
#   install.packages("Seurat")
# }
# library(Seurat)

# if (!"dplyr" %in% installed.packages()) {
#   install.packages("dplyr")
# }
# library(dplyr)

# if (!"patchwork" %in% installed.packages()) {
#   install.packages("patchwork")
# }
# library(patchwork)

# if (!"presto" %in% installed.packages()) {
#   install.packages('devtools')
#   devtools::install_github('immunogenomics/presto')
# }
# library(presto)

# # 增加 future.globals.maxSize 选项的值
# options(future.globals.maxSize = 800 * 1024^2)  # 800 MiB

# # 加载数据
# pbmc <- readRDS("data/omFAP.rds")

# # 创建Seurat对象
# #pbmc <- CreateSeuratObject(counts = pbmc$counts, project = "pbmc", meta.data = pbmc$meta.data)

# # 设置默认Assay为"RNA"
# DefaultAssay(pbmc) <- "RNA"

# # 计算线粒体基因的百分比
# #pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# # 可视化QC
# dir.create("figures")
# png("figures/VlnPlot.png", width = 800, height = 800)
# VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# dev.off()
# png("figures/ScatterPlot.png", width = 800, height = 800)
# plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
# dev.off()

# # 标准化
# pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# # Feature筛选
# pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# #确认Features
# top10 <- head(VariableFeatures(pbmc), 10)

# # 可视化
# png("figures/VariableFeaturePlot.png", width = 800, height = 800)
# plot1 <- VariableFeaturePlot(pbmc)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2
# dev.off()

# # Scaling
# all.genes <- rownames(pbmc)
# pbmc <- ScaleData(pbmc, features = all.genes)

# #Linear dimensional reduction (PCA)
# pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# #可视化
# png("figures/PCAPlot1.png", width = 800, height = 800)
# VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
# dev.off()
# png("figures/PCAPlot2.png", width = 800, height = 800)
# DimPlot(pbmc, reduction = "pca") + NoLegend()
# dev.off()
# png("figures/PCAPlotHeatmap.png", width = 800, height = 800)
# DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
# dev.off()
# #ElbowPlot
# png("figures/ElbowPlot.png", width = 800, height = 800)
# ElbowPlot(pbmc)
# dev.off()

# # Cluster the cells
# pbmc <- FindNeighbors(pbmc, dims = 1:10)
# pbmc <- FindClusters(pbmc, resolution = 0.5)
# head(Idents(pbmc), 5)

# # non-linear dimensional reduction (UMAP)
# pbmc <- RunUMAP(pbmc, dims = 1:10)
# png("figures/UMAPPlot.png", width = 800, height = 800)
# DimPlot(pbmc, reduction = "umap")
# dev.off()

# # Finding differentially expressed genes
# # REALLY TAKES TIME
# cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
# head(cluster2.markers, 5)
# cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
# head(cluster5.markers, 5)
# pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
# pbmc.markers %>%
#   group_by(cluster) %>%
#   dplyr::filter(avg_log2FC > 1)
# # save.image("tmp.RData")
# cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
# png("figures/Cluster0Markers.png", width = 800, height = 800)
# VlnPlot(pbmc, features = c("CD3E", "CD3D", "CD3G"), pt.size = 0.1) + NoLegend()
# dev.off()

library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
rm(list = ls())
setwd("test")
dir.create("figures")

# 读入数据 FROM 10X WEBSITE
pbmc <- Read10X("data/filtered_feature_bc_matrix")
pbmc <- CreateSeuratObject(pbmc)
glimpse(pbmc)

# QC
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(pbmc@assays$RNA))
HB.genes <- rownames(pbmc@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
pbmc[["percent.HB"]] <- PercentageFeatureSet(pbmc, features = HB.genes)
head(pbmc@meta.data)
col.num <- length(levels(pbmc@active.ident))
violin <- VlnPlot(pbmc,
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB"),
                 cols = rainbow(col.num),
                 pt.size = 0.1,
                 ncol = 4) +
                 theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axix.ticks.x = element_blank())
ggsave("figures/vlnplot_before_QC.pdf", plot = violin, width = 12, height = 6)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearlplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow = 1, legend = "none")
ggsave("figures/pearplot_before_QC.pdf", plot = pearlplot, width = 12, height = 5)

minGene = 300
maxGene = 4000
percentMT = 15

pbmc <- subset(pbmc, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < percentMT)
col.num <- length(levels(pbmc@active.ident))
violin <- VlnPlot(pbmc,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB"),
                  cols = rainbow(col.num),
                  pt.size = 0.1,
                  ncol = 4) +
                  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave("figures/vlnplot_after_QC.pdf", plot = violin, width = 12, height = 6)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
saveRDS(pbmc, file = "pbmc.RDS")
