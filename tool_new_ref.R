# source("default.R")

library(parallel)
library(tidyverse)
library(Seurat)
library(mclust)
library(MASS)
library(scID)
library(biomod2)

print("set parameter")
# input cell
object <- "all"

celltypes <- c("celltype", "celltype_merge")
celltype <- celltypes[2]


# For reference based
merges <- c("6metamodules", "4_merge_metamodules", "4_merge_metamodules_3celltypes", "4_merge_metamodules_mes")
merge <- merges[3]
# pick which celltype to be analyse

all_celltypes <- c("AClike", "MESlike", "NPClike", "OPClike")
chosen_celltypes <- all_celltypes[c(1, 3, 4)]


library(scID)
library(Seurat)

# Ref
neftel.smt <- readRDS("./../output/smrt_mal")
neftel.smt <- subset(neftel.smt, idents = chosen_celltypes)
reference_gem <- as.matrix(neftel.smt@assays$norm@data)
reference_clusters <- as.factor(neftel.smt@meta.data[, celltype])
names(reference_clusters) <- rownames(neftel.smt@meta.data)

# Target
gbm <- readRDS("./../output/seurat_gbm_qc")
gbm.list <- SplitObject(gbm, split.by = "split")
gbm.list[["run2_radiated_E31N"]] <- merge(gbm.list[["run2_radiated_E31N"]], y = gbm.list[["run1_radiated_E31N"]])
gbm.list[["run2_control_E31N"]] <- merge(gbm.list[["run2_control_E31N"]], y = gbm.list[["run1_control_E31N"]])
gbm.list[["run2_radiated_E26N"]] <- merge(gbm.list[["run2_radiated_E26N"]], y = gbm.list[["run1_radiated_E26N"]])
gbm.list[["run2_radiated_E24N"]] <- merge(gbm.list[["run2_radiated_E24N"]], y = gbm.list[["run1_radiated_E24N"]])
gbm.list[c("run1_radiated_E24N", "run1_radiated_E26N", "run1_control_E31N", "run1_radiated_E31N")] <- NULL

rm(gbm, neftel.smt)
gc()

scid_multiclass_edited <- function(target_gem = NULL, reference_gem = NULL,
                                   reference_clusters = NULL, markers = NULL,
                                   logFC = 0.5, normalize_reference = TRUE,
                                   estimate_weights_from_target = FALSE,
                                   weights = NULL, only_pos = FALSE) {
  # state the variable that will be use in the function
  markers <- markers.glob

  # ----------------------------------------------------------------------------------------------------
  # Data pre-processing
  # ----------------------------------------------------------------------------------------------------
  if (is.null(reference_gem) && is.null(reference_clusters) && is.null(markers)) {
    stop("Please provide either clustered reference data or list of markers for each reference cluster")
  }
  if (!is.null(reference_gem) && !is.null(reference_clusters)) {
    # Check all reference cells have a cluster ID
    common_cells <- intersect(names(reference_clusters), colnames(reference_gem))
    if (length(common_cells) == 0) {
      stop("None  of the reference cells has a cluster ID. Please check the reference_clusters list provided.")
    } else {
      reference_gem <- reference_gem[, common_cells]
      rownames(reference_gem) <- make.names(toupper(rownames(reference_gem)), unique = TRUE)

      # Remove genes that are zero across all cells
      reference_gem <- reference_gem[which(rowSums(reference_gem) != 0), ]
      reference_clusters <- reference_clusters[common_cells]
    }
  }

  if (!is.null(markers)) {
    # Check markers have gene and cluster columns
    if (length(intersect(c("gene", "cluster"), colnames(markers))) != 2) {
      stop("Please provide a data frame of markers with gene and cluster in columns")
    }
    markers$gene <- toupper(markers$gene)
  }


  # Target
  rownames(target_gem) <- make.names(toupper(rownames(target_gem)), unique = TRUE)
  # Remove genes that are zero across all cells
  target_gem <- target_gem[which(rowSums(target_gem) != 0), ]

  # ----------------------------------------------------------------------------------------------------
  # Stage 1: Find signature genes from reference data
  # ----------------------------------------------------------------------------------------------------
  if (is.null(markers)) {
    # Filter out signature genes that are not present in the target data
    markers <- markers[which(markers$gene %in% rownames(target_gem)), ]

    celltypes <- unique(markers$cluster)

    if (estimate_weights_from_target) {
      rm(reference_gem, reference_clusters)
    }
  } else {
    markers <- markers[which(markers$gene %in% rownames(target_gem)), ]
    celltypes <- unique(markers$cluster)
  }

  # Min-max normalization of target gem
  target_gem_norm <- t(apply(target_gem[unique(markers$gene), ], 1, function(x) normalize_gene(x)))
  target_gem_norm <- target_gem_norm[complete.cases(target_gem_norm), ]

  # ----------------------------------------------------------------------------------------------------
  # Stage 2: Weight signature genes
  # ----------------------------------------------------------------------------------------------------
  if (is.null(weights)) {
    if (estimate_weights_from_target) {
      weights <- list()
      for (i in 1:length(celltypes)) {
        celltype_markers <- markers[which(markers$cluster == celltypes[i]), ]
        positive_markers <- celltype_markers$gene[which(celltype_markers$avg_log2FC > 0)]
        negative_markers <- celltype_markers$gene[which(celltype_markers$avg_log2FC < 0)]
        training_groups <- choose_training_set(target_gem, positive_markers, negative_markers)
        signature_genes <- c(positive_markers, negative_markers)
        gene.weights <- scID_weight(target_gem_norm[signature_genes, , drop = FALSE], training_groups$in_pop, training_groups$out_pop)
        # If only positive markers are selected, truncate all negative weights to 0
        if (only_pos) {
          gene.weights[which(gene.weights < 0)] <- 0
        }
        # Make Inf weights 0
        gene.weights[is.infinite(gene.weights)] <- 0
        weights[[as.character(celltypes[i])]] <- gene.weights
      }
      names(weights) <- celltypes
    } else {
      if (!is.null(reference_gem) && !is.null(reference_clusters)) {
        weights <- list()
        # Normalize reference gem
        ref_gem_norm <- t(apply(reference_gem[unique(markers$gene), ], 1, function(x) normalize_gene(x)))
        ref_gem_norm <- ref_gem_norm[complete.cases(ref_gem_norm), ]
        for (i in 1:length(celltypes)) {
          signature_genes <- markers$gene[which(markers$cluster == celltypes[i])]
          true_cells <- names(reference_clusters)[which(reference_clusters == as.character(celltypes[i]))]
          false_cells <- setdiff(names(reference_clusters), true_cells)
          gene.weights <- scID_weight(gem = ref_gem_norm[signature_genes, , drop = FALSE], true_cells, false_cells)

          weights[[as.character(celltypes[i])]] <- gene.weights
          # If only positive markers are selected, truncate all negative weights to 0
          if (only_pos) {
            gene.weights[which(gene.weights < 0)] <- 0
          }
        }
        # Won't need reference data any more, remove for efficiency
        rm(reference_gem, reference_clusters, ref_gem_norm)
      } else {
        stop("Please provide reference data in order to calculate weights, choose to estimate weights from target data, or provide precompted gene weights.")
      }
    }
  }

  #----------------------------------------------------------------------------------------------------
  # Stage 3: Find scores and putative matches
  # ----------------------------------------------------------------------------------------------------

  scores <- data.frame(matrix(NA, length(celltypes), ncol(target_gem)), row.names = celltypes)
  colnames(scores) <- colnames(target_gem)

  full_scores <- data.frame(matrix(NA, length(celltypes), ncol(target_gem)), row.names = celltypes)
  colnames(full_scores) <- colnames(target_gem)

  for (i in 1:length(celltypes)) {
    celltype <- as.character(celltypes[i])
    signature <- intersect(names(weights[[celltype]]), rownames(target_gem_norm))
    weighted_gem <- weights[[celltype]][signature] * target_gem_norm[signature, , drop = FALSE]
    # Check if whole weighted gem is 0 (when all gene weighst are zero)
    if (all(weighted_gem == 0)) {
      full_scores[as.character(celltype), ] <- rep(0, ncol(full_scores))
      print(paste("oh no", i))
    } else {
      score <- colSums(weighted_gem) / sqrt(sum(weights[[celltype]]^2))

      matches <- final_populations_edited(score)

      scores[as.character(celltype), matches] <- scale(score[matches])
      full_scores[as.character(celltype), ] <- score
    }
  }

  # Resolve multi-class assignments

  labels <- apply(scores, 2, function(x) {
    ifelse(all(is.na(x)), "unassigned", rownames(scores)[which(x == max(x, na.rm = T))])
  })

  # Add convert unassigned to random assigned one
  # Check first whether we have unassigned or not
  unassigned_index <- which(labels == "unassigned")

  labels_edited <- labels
  full_scores_unassigned <- NA

  if (length(unassigned_index) > 0) {
    full_scores_unassigned <- full_scores[, unassigned_index, drop = FALSE]

    labels_unassigned <- apply(full_scores_unassigned, 2, function(x) {
      rownames(full_scores_unassigned)[which(x == max(x, na.rm = T))]
    })
    labels_edited[unassigned_index] <- labels_unassigned
  } else {
    message("no unassigned")
  }

  # return result
  list(scores = scores, full_scores = full_scores, labels = labels, labels_edited = labels_edited)
}

final_populations_edited <- function(score) {
  fit <- suppressMessages(mclust::densityMclust(score, G = 1:4, plot = FALSE, verbose = FALSE))

  # Calculate average scID score per group of cells
  avgScore <- rep(NA, length(unique(fit$classification)))
  names(avgScore) <- unique(fit$classification)
  for (ID in names(avgScore)) avgScore[ID] <- mean(score[names(which(fit$classification == ID))])

  matches <- names(fit$classification)[which(fit$classification == names(which(avgScore == max(avgScore))))]

  matches
}

find_markers_edited <- function(
    reference_gem, reference_clusters, logFC, only.pos,
    normalize_reference) {
  library(Seurat)
  so_ref <- CreateSeuratObject(reference_gem)
  if (normalize_reference) {
    so_ref <- suppressMessages(NormalizeData(so_ref, verbose = FALSE))
  }
  so_ref <- suppressMessages(ScaleData(so_ref, verbose = FALSE))
  Idents(so_ref) <- as.factor(reference_clusters)
  markers <- suppressMessages(FindAllMarkers(so_ref,
    test.use = "MAST",
    only.pos = only.pos, logfc.threshold = logFC, verbose = FALSE
  ))
  markers
}

markers.glob <- find_markers_edited(reference_gem, reference_clusters, logFC = 0.5, only.pos = FALSE, normalize_reference = FALSE)

# saveRDS(markers.glob,"./../output/markers_glob_scID")
# markers.glob <- readRDS("./../output/markers_glob_scID")

library(parallel)
# n.cores <- 2
n.cores <- parallel::detectCores()
# n.cores <- parallel::detectCores() - 3
all.tasks <- length(gbm.list)
i <- 0

# gbm.list.res <- lapply(X = gbm.list, FUN = function(x)  {


gbm.list.res <- mclapply(X = gbm.list, mc.cores = n.cores, FUN = function(x) {
  message("Running ", unique(x$split))
  i <<- i + 1
  message("it is ", i, " out of ", all.tasks)

  x <- NormalizeData(x, normalization.method = "RC", verbose = FALSE)
  x <- as.matrix(x@assays$RNA@data)
  # If want to set estimate_weights_from_target = TRUE, need older version of biomod2 (3.5.1)
  # scID_output <- tryCatch(scid_multiclass_edited(target_gem = x, reference_gem = reference_gem, reference_clusters = reference_clusters, logFC = 0.6, only_pos = FALSE,  estimate_weights_from_target = TRUE, normalize_reference = FALSE), error=function(e) NA)
  scID_output <- scid_multiclass_edited(target_gem = x, reference_gem = reference_gem, reference_clusters = reference_clusters, logFC = 0.6, only_pos = FALSE, estimate_weights_from_target = TRUE, normalize_reference = FALSE)

  #  gc()
  return(scID_output)
})

saveRDS(gbm.list.res, file = paste0("./../output/scID_", object, "_", merge, "_allassigned"))

res.df <- data.frame()
for (i in gbm.list.res) {
  res.each.df <- cbind(as.data.frame(i$labels, stringsAsFactors = TRUE), as.data.frame(i$labels_edited, stringsAsFactors = TRUE), data.frame(t(i$scores)))
  res.df <- rbind(res.df, res.each.df)
}
colnames(res.df)[1] <- "scID"
colnames(res.df)[2] <- "scID_edited"

# gbm.meta <- read.csv("./../output/gbm_meta.csv",row.names = 1)


# scid <- merge(res.df, gbm.meta, by = 'row.names', all = TRUE)

write.csv(res.df, paste0("./../output/scID_", object, "_", merge, "_allassigned", ".csv"), row.names = TRUE)
