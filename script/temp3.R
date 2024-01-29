function (target_gem = NULL, reference_gem = NULL, reference_clusters = NULL, 
          markers = NULL, logFC = 0.5, normalize_reference = TRUE, 
          estimate_weights_from_target = FALSE, weights = NULL, only_pos = FALSE) 
{
  if (is.null(reference_gem) && is.null(reference_clusters) && 
      is.null(markers)) {
    stop("Please provide either clustered reference data or list of markers for each reference cluster")
  }
  if (!is.null(reference_gem) && !is.null(reference_clusters)) {
    common_cells <- intersect(names(reference_clusters), 
                              colnames(reference_gem))
    if (length(common_cells) == 0) {
      stop("None  of the reference cells has a cluster ID. Please check the reference_clusters list provided.")
    }
    else {
      reference_gem <- reference_gem[, common_cells]
      rownames(reference_gem) <- make.names(toupper(rownames(reference_gem)), 
                                            unique = TRUE)
      reference_gem <- reference_gem[which(rowSums(reference_gem) != 
                                             0), ]
      reference_clusters <- reference_clusters[common_cells]
    }
  }
  if (!is.null(markers)) {
    if (length(intersect(c("gene", "cluster"), colnames(markers))) != 
        2) {
      stop("Please provide a data frame of markers with gene and cluster in columns")
    }
    markers$gene <- toupper(markers$gene)
  }
  rownames(target_gem) <- make.names(toupper(rownames(target_gem)), 
                                     unique = TRUE)
  target_gem <- target_gem[which(rowSums(target_gem) != 0), 
  ]
  if (is.null(markers)) {
    message("Stage 1: extract signatures genes from reference clusters")
    markers <- find_markers(reference_gem, reference_clusters, 
                            logFC, only.pos = only_pos, normalize_reference = normalize_reference)
    markers <- markers[which(markers$gene %in% rownames(target_gem)), 
    ]
    celltypes <- unique(markers$cluster)
    if (estimate_weights_from_target) {
      rm(reference_gem, reference_clusters)
    }
  }
  else {
    markers <- markers[which(markers$gene %in% rownames(target_gem)), 
    ]
    celltypes <- unique(markers$cluster)
  }
  target_gem_norm <- t(apply(target_gem[unique(markers$gene), ], 1, function(x) normalize_gene(x)))
  target_gem_norm <- target_gem_norm[complete.cases(target_gem_norm), 
  ]
  if (is.null(weights)) {
    if (estimate_weights_from_target) {
      message("Stage 2: Estimate weights of signature genes from target")
      weights <- list()
      for (i in 1:length(celltypes)) {
        celltype_markers <- markers[which(markers$cluster == 
                                            celltypes[i]), ]
        positive_markers <- celltype_markers$gene[which(celltype_markers$avg_log2FC > 
                                                          0)]
        negative_markers <- celltype_markers$gene[which(celltype_markers$avg_log2FC < 
                                                          0)]
        training_groups <- choose_training_set(target_gem, 
                                               positive_markers, negative_markers)
        signature_genes <- c(positive_markers, negative_markers)
        gene.weights <- scID_weight(target_gem_norm[signature_genes, 
                                                    , drop = FALSE], training_groups$in_pop, training_groups$out_pop)
        if (only_pos) {
          gene.weights[which(gene.weights < 0)] <- 0
        }
        gene.weights[is.infinite(gene.weights)] <- 0
        weights[[as.character(celltypes[i])]] <- gene.weights
        svMisc::progress(i * 100/length(celltypes))
        Sys.sleep(0.01)
        if (i == length(celltypes)) 
          cat("Done!")
      }
      names(weights) <- celltypes
    }
    else {
      if (!is.null(reference_gem) && !is.null(reference_clusters)) {
        message("Stage 2: Estimate weights of signature genes from reference")
        weights <- list()
        ref_gem_norm <- t(apply(reference_gem[unique(markers$gene), 
        ], 1, function(x) normalize_gene(x)))
        ref_gem_norm <- ref_gem_norm[complete.cases(ref_gem_norm), 
        ]
        for (i in 1:length(celltypes)) {
          signature_genes <- markers$gene[which(markers$cluster == 
                                                  celltypes[i])]
          true_cells <- names(reference_clusters)[which(reference_clusters == 
                                                          as.character(celltypes[i]))]
          false_cells <- setdiff(names(reference_clusters), 
                                 true_cells)
          gene.weights <- scID_weight(gem = ref_gem_norm[signature_genes, 
                                                         , drop = FALSE], true_cells, false_cells)
          weights[[as.character(celltypes[i])]] <- gene.weights
          if (only_pos) {
            gene.weights[which(gene.weights < 0)] <- 0
          }
          svMisc::progress(i * 100/length(celltypes))
          Sys.sleep(0.01)
          if (i == length(celltypes)) 
            cat("Done!")
        }
        rm(reference_gem, reference_clusters, ref_gem_norm)
      }
      else {
        stop("Please provide reference data in order to calculate weights, choose to estimate weights from target data, or provide precompted gene weights.")
      }
    }
  }
  message("Stage 3.1-2: Calculate scores and find matching cells")
  scores <- data.frame(matrix(NA, length(celltypes), ncol(target_gem)), 
                       row.names = celltypes)
  colnames(scores) <- colnames(target_gem)
  full_scores <- data.frame(matrix(NA, length(celltypes), ncol(target_gem)), 
                            row.names = celltypes)
  colnames(full_scores) <- colnames(target_gem)
  for (i in 1:length(celltypes)) {
    celltype <- as.character(celltypes[i])
    signature <- intersect(names(weights[[celltype]]), rownames(target_gem_norm))
    weighted_gem <- weights[[celltype]][signature] * target_gem_norm[signature, 
                                                                     , drop = FALSE]
    if (all(weighted_gem == 0)) {
      full_scores[as.character(celltype), ] <- rep(0, ncol(full_scores))
    }
    else {
      score <- colSums(weighted_gem)/sqrt(sum(weights[[celltype]]^2))
      matches <- final_populations(score)
      scores[as.character(celltype), matches] <- scale(score[matches])
      full_scores[as.character(celltype), ] <- score
    }
    if (i == length(celltypes)) 
      cat("Done!")
  }
  message("Stage 3.3: Resolve multiclass assignments")
  labels <- apply(scores, 2, function(x) {
    ifelse(all(is.na(x)), "unassigned", rownames(scores)[which(x == 
                                                                 max(x, na.rm = T))])
  })
  list(markers = markers, estimated_weights = weights, labels = labels, 
       scores = full_scores)
}

function (score) 
{
  fit <- mclust::densityMclust(score)
  avgScore <- rep(NA, length(unique(fit$classification)))
  names(avgScore) <- unique(fit$classification)
  for (ID in names(avgScore)) avgScore[ID] <- mean(score[names(which(fit$classification == 
                                                                       ID))])
  matches <- names(fit$classification)[which(fit$classification == 
                                               names(which(avgScore == max(avgScore))))]
  matches
}


function (reference_gem, reference_clusters, logFC, only.pos, 
          normalize_reference) 
{
  library(Seurat)
  so_ref <- CreateSeuratObject(reference_gem)
  if (normalize_reference) {
    so_ref <- suppressMessages(NormalizeData(so_ref))
  }
  so_ref <- suppressMessages(ScaleData(so_ref))
  Idents(so_ref) <- as.factor(reference_clusters)
  markers <- suppressMessages(FindAllMarkers(so_ref, test.use = "MAST", 
                                             only.pos = only.pos, logfc.threshold = logFC))
  markers
}