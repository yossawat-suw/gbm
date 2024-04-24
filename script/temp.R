#processing and visualization
.libPaths("/home/point/R/x86_64-pc-linux-gnu-library/4.3")
setwd("~/Documents/work/research/pmayp/Project/lucy/gbm")

library(pacman)
p_load(qs,Seurat)
seurat_obj <- qread( "output/seurat_objects/neurodevelopmental_filtered.qs")



seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(seurat_obj)



gc()




seurat_obj <- ScaleData(seurat_obj, features = all.genes)



gc()



seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)




# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(seurat_obj), 10)
#
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(seurat_obj)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2




gc()



seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))



gc()



seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

.libPaths()

gc()

qsave(seurat_obj, "output/seurat_objects/neurodevelopmental_filtered_processed.qs")
