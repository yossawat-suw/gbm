#library(escape)
library(Seurat)

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



gbm.list_positive_TRUE <- readRDS(paste0("output/deg_gsea/gbm_escaped_normpos/",args,".rds"))

print("run_find_marker")


Idents(gbm.list_positive_TRUE) <- "radiation"
all.markers <- FindAllMarkers(gbm.list_positive_TRUE,
               assay = "escape.ssGSEA_normalized",
               min.pct = 0,
               logfc.threshold = 0,
               verbose = TRUE)



saveRDS(all.markers,paste0("output/deg_gsea/all_marker/",args,".rds"))




