#library(escape)
#library(Seurat)

# .libPaths("/home/point/Documents/work/research/pmayp/Project/lucy/gbm/renv/library/R-4.1/x86_64-pc-linux-gnu")
# .libPaths("/home/point/R/x86_64-pc-linux-gnu-library/4.3")


#commandArgs picks up the variables you pass from the command line
args <- commandArgs(trailingOnly = TRUE)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("escape")

setwd(here::here())


source("script/my_escape.R")


print("Load GS")

gs_all <- readRDS(file = "output/deg_gsea/gs_all.rds")

print("Load GBM")
gbm.list_positive_TRUE <- readRDS(paste0("output/deg_gsea/gbm_escaped/",args,".rds"))


print("Differential GSEA")

print("run_norm")
gbm.list_positive_TRUE <-  performNormalization_edited(gbm.list_positive_TRUE, 
                              assay = "escape.ssGSEA",
                              gene.sets = gs_all,make.positive = FALSE)


saveRDS(gbm.list_positive_TRUE,paste0("output/deg_gsea/gbm_escaped_normpos/",args,".rds"))

# print("run_find_marker")
# gc()
# 
# #plan(multicore,workers = parallel::detectCores() -3)
# plan(multicore,workers = 2)
# all.markers.list <-  future_lapply(gbm.list_positive_TRUE,function(x) {
#   Idents(x) <- "radiation"
#   FindAllMarkers(x, 
#                  assay = "escape.ssGSEA_normalized", 
#                  min.pct = 0,
#                  logfc.threshold = 0)
# },future.seed = 7)
# 
# 
# saveRDS(all.markers.list,"output/deg_gsea/gbm_gsea_findallmarker.rds")




