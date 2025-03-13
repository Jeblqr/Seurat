# 安装并导入必要的包
if (!"Seurat" %in% installed.packages()) {
  install.packages("Seurat")
}
library(Seurat)

if (!"dplyr" %in% installed.packages()) {
  install.packages("dplyr")
}
library(dplyr)

if (!"patchwork" %in% installed.packages()) {
  install.packages("patchwork")
}
library(patchwork)

# 增加 future.globals.maxSize 选项的值
options(future.globals.maxSize = 800 * 1024^2)  # 800 MiB

# 加载数据
omFAP <- readRDS("data/omFAP.rds")

# 创建Seurat对象
#omFAP <- CreateSeuratObject(counts = omFAP$counts, project = "omFAP", meta.data = omFAP$meta.data)

# 设置默认Assay为"RNA"
DefaultAssay(omFAP) <- "RNA"

# 计算线粒体基因的百分比
#omFAP[["percent.mt"]] <- PercentageFeatureSet(omFAP, pattern = "^MT-")

# 可视化QC
dir.create("figures")
png("figures/VlnPlot.png", width = 800, height = 800)
VlnPlot(omFAP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
png("figures/ScatterPlot.png", width = 800, height = 800)
plot1 <- FeatureScatter(omFAP, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(omFAP, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

# 标准化
omFAP <- NormalizeData(omFAP, normalization.method = "LogNormalize", scale.factor = 10000)

# Feature筛选
omFAP <- FindVariableFeatures(omFAP, selection.method = "vst", nfeatures = 2000)

#确认Features
top10 <- head(VariableFeatures(omFAP), 10)

# 可视化
png("figures/VariableFeaturePlot.png", width = 800, height = 800)
plot1 <- VariableFeaturePlot(omFAP)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

# Scaling
all.genes <- rownames(omFAP)
omFAP <- ScaleData(omFAP, features = all.genes)

#PCA
omFAP <- RunPCA(omFAP, features = VariableFeatures(object = omFAP))

#可视化
png("figures/PCAPlot1.png", width = 800, height = 800)
VizDimLoadings(omFAP, dims = 1:2, reduction = "pca")
dev.off()
png("figures/PCAPlot2.png", width = 800, height = 800)
DimPlot(omFAP, reduction = "pca") + NoLegend()
dev.off()
png("figures/PCAPlotHeatmap.png", width = 800, height = 800)
DimHeatmap(omFAP, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

#ElbowPlot
png("figures/ElbowPlot.png", width = 800, height = 800)
ElbowPlot(omFAP)
dev.off()