.libPaths("/home/point/R/x86_64-pc-linux-gnu-library/4.3")

library(escape)
library(Seurat)
library(parallel)
library(future)
library(future.apply)
library("matrixStats")
library(pbapply)
library(BiocParallel)

setwd(here::here())


source("script/my_escape.R")

# print("prepare gs")
# #https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
# GS.msigdb <- getGeneSets(library = c("H","C2","C3","C4","C5","C6"))
# 
# data("escape.gene.sets", package="escape")
# gene.sets <- escape.gene.sets
# sen_lucy <- read.csv("data/lucy_senesence_genes.csv")
# sen_lucy <- colnames(sen_lucy)
# 
# sen_mayo <- readxl::read_xlsx("data/gene_set/senmayo.xlsx")
# sen_mayo <- sen_mayo$`Gene(human)`
# gen.sets_manual <- list(sen_lucy = sen_lucy,
#                         sen_mayo = sen_mayo)
# gs_all <- c(GS.msigdb,gene.sets,gen.sets_manual)
# 
# names(gs_all) <- stringr::str_replace_all(names(gs_all), "_", "-")
# 
# saveRDS(gs_all,file = "output/deg_gsea/gs_all.rds")
# gc()

gs_all <- readRDS(file = "output/deg_gsea/gs_all.rds")
print("prepare GBM")

# gbm <- readRDS("output/seurat_objects/seurat_gbm_qc")
# all_genes <- rownames(gbm)
# 
# gbm.list <- SplitObject(gbm, split.by = "donor_id")


print("runEscape")
# gbm.list <- pblapply(gbm.list, function(x) {
#   x <-runEscape(x,
#                 method = "ssGSEA",
#                 gene.sets = gs_all,
#                 groups = 1000,
#                 min.size = 5,
#                 new.assay.name = "escape.ssGSEA",
#                 BPPARAM = MulticoreParam(workers = (parallel::detectCores() - 3))
#   )
#   return(x)
# })

#plan(multicore,workers = parallel::detectCores() -3) too expensive
# plan(multicore,workers = 4)
# 
# gbm.list <- future_lapply(gbm.list, function(x) {
#   x <-runEscape(x,
#                 method = "ssGSEA",
#                 gene.sets = gs_all,
#                 groups = 1000,
#                 min.size = 5,
#                 new.assay.name = "escape.ssGSEA"
#   )
#   return(x)
# },future.seed = 7)
# 
# 
# gc()
# print("save Escape")
# saveRDS(gbm.list,"output/deg_gsea/gbm_seurat_list_escaped.rds")

# print("Load GBM")
# 
# gbm.list <- readRDS("output/deg_gsea/gbm_seurat_list_escaped.rds")

# gbm.list <- future_lapply(gbm.list, function(x) {
#   performNormalization_edited(x,
#                               assay = "escape.ssGSEA",
#                               gene.sets = gs_all)
# },future.seed = 7)
# gc()
# 
# 
# 
# 

# 
# gbm.combined <- Reduce(merge,gbm.list)
# 
# saveRDS(gbm.combined,file = "output/deg_gsea/gbm_seurat_processed_normalized_combined.rds")


gbm.list <- readRDS("output/deg_gsea/gbm_seurat_list_escaped.rds")
for (i in seq_along(gbm.list)) {
  saveRDS(gbm.list[[i]],paste0("output/deg_gsea/gbm_escaped/"),names(gbm.list)[i])
}

print("Differential GSEA")
gc()
#Differential GSEA
#plan(multicore,workers = parallel::detectCores() -3)
plan(multicore,workers = 2)
print("run_norm")
gbm.list_positive_TRUE <-  future_lapply(gbm.list,function(x) {
  performNormalization_edited(x, 
                              assay = "escape.ssGSEA",
                              gene.sets = gs_all,make.positive = FALSE)
},future.seed = 7)

saveRDS(gbm.list_positive_TRUE,"output/deg_gsea/gbm_seurat_list_escaped_normalized_positive_true.rds")

print("run_find_marker")
gc()

#plan(multicore,workers = parallel::detectCores() -3)
plan(multicore,workers = 2)
all.markers.list <-  future_lapply(gbm.list_positive_TRUE,function(x) {
  Idents(x) <- "radiation"
  FindAllMarkers(x, 
                 assay = "escape.ssGSEA_normalized", 
                 min.pct = 0,
                 logfc.threshold = 0)
},future.seed = 7)


saveRDS(all.markers.list,"output/deg_gsea/gbm_gsea_findallmarker.rds")




