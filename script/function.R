sctype_score_edited <- function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...) {
  # check input matrix
  if (!is.matrix(scRNAseqData)) {
    warning("scRNAseqData doesn't seem to be a matrix")
  } else {
    if (sum(dim(scRNAseqData)) == 0) {
      warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
    }
  }

  # marker sensitivity
  marker_stat <- sort(table(unlist(gs)), decreasing = T)
  marker_sensitivity <- data.frame(
    score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0, 1), from = c(length(gs), 1)),
    gene_ = names(marker_stat), stringsAsFactors = !1
  )

  # convert gene names to Uppercase
  if (gene_names_to_uppercase) {
    rownames(scRNAseqData) <- toupper(rownames(scRNAseqData))
  }

  # subselect genes only found in data
  names_gs_cp <- names(gs)
  names_gs_2_cp <- names(gs2)
  gs <- lapply(1:length(gs), function(d_) {
    GeneIndToKeep <- rownames(scRNAseqData) %in% as.character(gs[[d_]])
    rownames(scRNAseqData)[GeneIndToKeep]
  })
  gs2 <- lapply(1:length(gs2), function(d_) {
    GeneIndToKeep <- rownames(scRNAseqData) %in% as.character(gs2[[d_]])
    rownames(scRNAseqData)[GeneIndToKeep]
  })
  names(gs) <- names_gs_cp
  names(gs2) <- names_gs_2_cp
  cell_markers_genes_score <- marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)), ]

  # z-scale if not
  if (!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData

  # multiple by marker sensitivity
  for (jj in 1:nrow(cell_markers_genes_score)) {
    Z[cell_markers_genes_score[jj, "gene_"], ] <- Z[cell_markers_genes_score[jj, "gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }

  # subselect only with marker genes
  Z <- Z[unique(c(unlist(gs), unlist(gs2))), ]

  # combine scores
  es <- do.call("rbind", lapply(names(gs), function(gss_) {
    sapply(1:ncol(Z), function(j) {
      gs_z <- Z[gs[[gss_]], j]
      gz_2 <- Z[gs2[[gss_]], j] * -1
      sum_t1 <- (sum(gs_z) / sqrt(length(gs_z)))
      sum_t2 <- sum(gz_2) / sqrt(length(gz_2))
      if (is.na(sum_t2)) {
        sum_t2 <- 0
      }
      sum_t1 + sum_t2
    })
  }))

  dimnames(es) <- list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all), , drop = FALSE] # remove na rows

  es.max
}



adjust_composition_center <- function(vec, adjustment_factor) {
  if(!is.numeric(vec)) {
    stop("Input must be a numeric vector of length 3")
  }
  
  # Calculate the average of the vector
  average_value <- mean(vec)
  
  # Adjust each value towards the average
  vec <- vec + adjustment_factor * (average_value - vec)
  
  # Normalize the vector so that the sum equals 1
  vec <- vec / sum(vec)
  
  return(vec)
}


scid_multiclass_edited <- function(target_gem = NULL, reference_gem = NULL,
                                   reference_clusters = NULL, markers = NULL,
                                   logFC = 0.5, normalize_reference = TRUE,
                                   estimate_weights_from_target = FALSE,
                                   weights = NULL, only_pos = FALSE,
                                   markers.glob = markers.glob) {
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
      max_value <- max(x, na.rm = TRUE)
      tie_indices <- which(x == max_value)
      if (length(tie_indices) > 1) {
        "unassigned"  # Assign "unassigned" when tie condition occurs
      } else {
        rownames(full_scores_unassigned)[tie_indices]
      }
    })
    
    # check whether the labels_unassigned and the same name/order of labels_edited[unassigned_index]
    if (sum(names(labels_edited[unassigned_index]) == names(labels_unassigned)) == length(labels_unassigned)) {
      labels_edited[unassigned_index] <- unname(labels_unassigned)
    } else {
      message("not equal ")
    }
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

#the edited 
find_markers_edited <- function(
    reference_gem, reference_clusters, logFC, only.pos,
    normalize_reference) {
  library(Seurat)
  so_ref <- CreateSeuratObject(reference_gem) 
  so_ref@assays$RNA$data <- so_ref@assays$RNA$counts ## <- edit here as the count is already norm data
  so_ref <- DietSeurat(so_ref,layers = c("data"))
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

find_markers_edited_norm <- function(
    reference_gem, reference_clusters, logFC, only.pos,
    normalize_reference) {
  library(Seurat)
  so_ref <- CreateSeuratObject(reference_gem) 
  if (normalize_reference) {
    so_ref <- suppressMessages(NormalizeData(so_ref, normalization.method = "RC", verbose = FALSE,scale.factor = 1e6))
  } 
  #so_ref <- suppressMessages(ScaleData(so_ref, verbose = FALSE))
  Idents(so_ref) <- as.factor(reference_clusters)
  markers <- suppressMessages(FindAllMarkers(so_ref,
                                             test.use = "MAST",
                                             only.pos = only.pos, logfc.threshold = logFC, verbose = FALSE
  ))
  markers
}
##original find_marker 
# function (reference_gem, reference_clusters, logFC, only.pos, 
#           normalize_reference) 
# {
#   library(Seurat)
#   so_ref <- CreateSeuratObject(reference_gem)
#   if (normalize_reference) {
#     so_ref <- suppressMessages(NormalizeData(so_ref))
#   }
#   so_ref <- suppressMessages(ScaleData(so_ref))
#   Idents(so_ref) <- as.factor(reference_clusters)
#   markers <- suppressMessages(FindAllMarkers(so_ref, test.use = "MAST", 
#                                              only.pos = only.pos, logfc.threshold = logFC))
#   markers
# }

#For Scsorter

#' Preprocess Data
data_preprocess_edited = function(expr, anno_processed) {
  designmat = anno_processed[[1]]
  weightmat = anno_processed[[2]]
  
  rownames(expr) = toupper(rownames(expr))
  rownames(designmat) = toupper(rownames(designmat))
  
  markers = rownames(designmat)
  
  markers_avail = markers %in% rownames(expr)
  
  if (sum(markers_avail) < length(markers)) {
    warning(paste('The following specified marker genes are not found from the expression data: ',
                  paste(markers[!markers_avail], collapse = ', '), '.', sep = ''))
    
    designmat = designmat[markers_avail, ]
    csdmat = colSums(designmat)
    
    if(sum(csdmat > 0) < ncol(designmat)) {
      stop(paste('None of the marker genes specfied for the following cell types could be found from the expression data: ',
                 paste(colnames(designmat)[csdmat == 0], collapse = ', '), '. \n Please maker sure the marker gene names are correctly specified.', sep = ''))
    }
  }
  
  dmat_gene_names = rownames(designmat)
  dat_gene_names = rownames(expr)
  
  picker = dat_gene_names %in% dmat_gene_names
  
  expr_mk = expr[picker, ]
  expr_rt = expr[!picker,]
  
  #reorder genes so that the order of expr mat and design mat matches
  ror = rep(0, nrow(designmat))
  for(i in 1:nrow(designmat)){
    ror[i] = which(rownames(expr_mk) == dmat_gene_names[i])
  }
  expr_mk = expr_mk[ror,]
  
  expr_cb = rbind(expr_mk, expr_rt)
  
  return(list(dat = expr_cb, designmat = designmat, weightmat = weightmat))
}

scSorter_edited <- function (expr, anno, default_weight = 2, n_start = 10, alpha = 0, 
          u = 0.05, max_iter = 100, setseed = 0) 
{
  anno_processed = design_matrix_builder(anno, default_weight)
  dt = data_preprocess_edited(expr, anno_processed)
  dat = dt[[1]]
  designmat = dt[[2]]
  weightmat = dt[[3]]
  c_cost = NULL
  c_mu = list()
  c_clus = list()
  for (i in 1:n_start) {
    set.seed(i + setseed)
    t1 = Sys.time()
    pred_ot = update_func(as.matrix(dat), designmat, weightmat, 
                          unknown_threshold1 = alpha, unknown_threshold2 = u, 
                          max_iter = max_iter)
    t2 = Sys.time()
    c_cost = c(c_cost, pred_ot[[3]])
    c_mu[[i]] = pred_ot[[1]]
    c_clus[[i]] = pred_ot[[2]]
  }
  pk = which.min(c_cost)
  pred_clus = c_clus[[pk]]
  pred_clus = c(colnames(designmat), rep("Unknown", ncol(designmat)))[pred_clus]
  pred_mu = c_mu[[pk]]
  return(list(Pred_Type = pred_clus, Pred_param = pred_mu))
}


calculate_kappa_for_combination <- function(data, tools) {
  # Subset the data for the given combination of tools
  subset_data <- data[, tools, drop = FALSE]
  
  # Convert to count format (adjust based on your data's structure)
  
  # Calculate and return Fleiss' Kappa
  kappa_result <- kappam.fleiss(subset_data)
  return(kappa_result$value)
}


# Row is sample id, col is celltype 
clr.modified <- function(x, base = exp(1), pseudo_count = 1e-5) {
  # Adding pseudo-count to zeros
  x[x == 0] <- pseudo_count
  
  if (dim(x)[2] == 1) {
    res <- list(x.clr = x, gm = rep(1, dim(x)[1]))
  } else {
    geometricmean <- function(x) {
      if (any(na.omit(x == 0))) {
        0
      } else {
        exp(mean(log(unclass(x)[is.finite(x) & x > 0])))
      }
    }
    gm <- apply(x, 1, geometricmean)
    x.clr <- log(x / gm, base)
    res <- x.clr
  }
  return(res)
}


geometric_mean <- function(x) {
  exp(mean(log(x)))
}



checkInputFormat_edited <- function (waspas_db) 
{
  if (missing(waspas_db)) 
    return("Parameter waspas_db is missing")
  if (!is.data.frame(waspas_db)) 
    return("Parameter waspas_db must be a data.frame")
  tryCatch({
    proc_step <- "Indicators"
    indicators <- c("C", "F", "W")
    i1 <- toupper(substr(toString(waspas_db[1, 1]), 1, 1))
    i2 <- toupper(substr(toString(waspas_db[2, 1]), 1, 1))
    i3 <- toupper(substr(toString(waspas_db[3, 1]), 1, 1))
    indicators_uploaded <- sort(c(i1, i2, i3))
    if (!identical(indicators, indicators_uploaded)) 
      stop()
    flags <- sliceData(waspas_db, "F")
    just_bc <- sort(unique(toupper(substr(flags, 1, 1))))
    proc_step <- "Flags"
    # if (!identical(just_bc, c("B", "C"))) 
    #     stop()
    if (!any(just_bc %in% c("B", "C")))
      stop("Invalid first characters in 'F' column. Only 'B' or 'C' are allowed.")
    
    proc_step <- "Weights-1"
    weights <- sliceData(waspas_db, "W")
    weights <- sapply(weights, as.numeric)
    proc_step <- "Weights-2"
    #if (sum(weights) != 1) 
    #stop()
    proc_step <- "Values"
    values <- sliceData(waspas_db, "V")
    values <- sapply(values, as.numeric)
    if (sum(is.na(values)) > 0) 
      stop()
    proc_step <- "End"
  }, error = function(cond) {
    stop(paste("E[CI]", cond))
  }, warning = function(cond) {
    if (grepl("NAs intro", cond)) {
    }
  }, finally = {
    if (proc_step == "Indicators") {
      return(paste("Error: Check the indicators in cells [1:3, 1], the strings", 
                   " must be intiated with 'C', 'F' or 'W'"))
    }
    else if (proc_step == "Flags") {
      return(paste("Error: Vector of flags must contains just strings initiated", 
                   "with B or C (i.e. b, c, B, C, Cost, Benefit, Ben etc.)"))
    }
    else if (proc_step == "Weights-1") {
      return("Error: Check Weights values, all must be numeric")
    }
    else if (proc_step == "Weights-2") {
      return("Error: Values in Vector of Weights must summarize 1")
    }
    else if (proc_step == "Values") {
      return("Error: Check Aternatives x Criteria values, all must be numeric")
    }
    else {
      return(TRUE)
    }
  })
}

waspasR_edited <- function (waspas_df, lambda) 
{
  if (missing(waspas_df)) 
    return("Parameter waspas_df is missing")
  if (missing(lambda)) 
    return("Parameter lambda is missing")
  format_ok <- checkInputFormat_edited(waspas_df)
  if (is.character(format_ok)) 
    return(format_ok)
  Alternative <- sliceData(waspas_df, "A")
  criteria <- sliceData(waspas_df, "C")
  weights <- sliceData(waspas_df, "W")
  flags <- sliceData(waspas_df, "F")
  values <- sliceData(waspas_df, "V")
  normalized <- normalize_edited(values, flags)
  wsm <- calcWSM_edited(normalized, weights)
  wpm <- calcWPM_edited(normalized, weights)
  waspas <- applyLambda(wsm, wpm, lambda)
  waspas_matrix <- data.frame(matrix(nrow = nrow(waspas_df) - 
                                       1, ncol = ncol(waspas_df) + 3))
  colnames(waspas_matrix) <- cbind("Alternative", criteria, 
                                   "WSM_Rank", "WPM_Rank", "WASPAS_Rank")
  qtd_rows <- nrow(waspas_matrix)
  qtd_cols_weights <- ncol(weights) + 1
  qtd_cols_flags <- ncol(flags) + 1
  qtd_cols_values <- ncol(values) + 1
  waspas_matrix[1, 1] <- "W"
  waspas_matrix[1, 2:qtd_cols_weights] <- weights
  waspas_matrix[2, 1] <- "F"
  waspas_matrix[2, 2:qtd_cols_flags] <- flags
  waspas_matrix[3:qtd_rows, 1] <- t(Alternative)
  waspas_matrix[3:qtd_rows, 2:qtd_cols_values] <- values
  waspas_matrix[3:qtd_rows, "WSM_Rank"] <- waspas[, "WSM_Rank"]
  waspas_matrix[3:qtd_rows, "WPM_Rank"] <- waspas[, "WPM_Rank"]
  waspas_matrix[3:qtd_rows, "WASPAS_Rank"] <- waspas[, "WASPAS_Rank"]
  return(as.data.frame(waspas_matrix))
}


normalize_edited <- function (normalized_matrix, vec_cost_benefit) 
{
  tryCatch({
    sapply(normalized_matrix, as.numeric)
    if (length(vec_cost_benefit) != ncol(normalized_matrix)) {
      return(paste("Error: The cost - benefit flags array must be the same size", 
                   "as the number of criteria"))
    }
    just_bc <- sort(unique(toupper(substr(vec_cost_benefit, 
                                          1, 1))))
    # if (!identical(just_bc, c("B", "C"))) {
    #     return(paste("Error: Vector of flags must contains just strings initiated", 
    #         "with B or C (i.e. b, c, B, C, Cost, Benefit, Ben etc.)"))
    # }
    if (!any(just_bc %in% c("B", "C")))
      stop("Invalid first characters in 'F' column. Only 'B' or 'C' are allowed.")
    total_rows <- nrow(normalized_matrix)
    flags <- toupper(substr(vec_cost_benefit, 1, 1))
    for (iCol in seq_len(ncol(normalized_matrix))) {
      alternative_vals <- normalized_matrix[1:total_rows, 
                                            iCol]
      alternative_vals <- sapply(alternative_vals, as.numeric)
      maxv <- max(alternative_vals)
      minv <- min(alternative_vals)
      for (iRow in seq_len(nrow(normalized_matrix))) {
        if (flags[iCol] == "C") {
          normalized_matrix[iRow, iCol] <- toString(minv/as.numeric(normalized_matrix[iRow, 
                                                                                      iCol]))
        }
        else {
          normalized_matrix[iRow, iCol] <- toString(as.numeric(normalized_matrix[iRow, 
                                                                                 iCol])/maxv)
        }
      }
    }
    return(as.data.frame(normalized_matrix))
  }, error = function(cond) {
    stop(paste("E[N]", cond))
  }, warning = function(cond) {
    if (grepl("NAs intro", cond)) {
      return("W[N] Error: Some non numeric-alike value was found")
    }
  })
}

calcWSM_edited <- function (normal_db, vec_weights) 
{
  tryCatch({
    if (length(vec_weights) != ncol(normal_db)) {
      return(paste("Error: The weight vector must be the same size as the", 
                   "number of criteria"))
    }
    # if (sum(sapply(vec_weights, as.numeric)) != 1) {
    #     return("Error: Values in Vector of Weights must summarize 1")
    # }
    WSM_Rank <- rep(0, nrow(normal_db))
    Alternative <- seq_len(nrow(normal_db))
    wsm <- cbind(Alternative, WSM_Rank)
    for (iCol in seq_len(ncol(normal_db))) {
      for (iRow in seq_len(nrow(normal_db))) {
        normal_db[iRow, iCol] <- toString(as.numeric(normal_db[iRow, 
                                                               iCol]) * as.numeric(vec_weights[iCol]))
      }
    }
    for (iRow in seq_len(nrow(normal_db))) {
      wsm[iRow, "WSM_Rank"] <- sum(sapply(normal_db[iRow, 
      ], as.numeric))
    }
    wsm_db <- wsm[, c("Alternative", "WSM_Rank")]
    return(wsm_db)
  }, error = function(cond) {
    stop(paste("E[S]", cond))
  }, warning = function(cond) {
    if (grepl("NAs intro", cond)) {
      return("W[S] Error: Some non numeric - alike value was found")
    }
  })
}

calcWPM_edited <- function (normal_db, vec_weights) 
{
  tryCatch({
    if (length(vec_weights) != ncol(normal_db)) {
      return(paste("Error: The weight vector must be the same size as", 
                   "the number of criteria"))
    }
    # if (sum(sapply(vec_weights, as.numeric)) != 1) {
    #     return("Error: Values in Vector of Weights must summarize 1")
    # }
    WPM_Rank <- rep(0, nrow(normal_db))
    Alternative <- seq_len(nrow(normal_db))
    wpm <- cbind(Alternative, WPM_Rank)
    for (iCol in seq_len(ncol(normal_db))) {
      for (iRow in seq_len(nrow(normal_db))) {
        normal_db[iRow, iCol] <- toString(as.numeric(normal_db[iRow, 
                                                               iCol])^as.numeric(vec_weights[iCol]))
      }
    }
    for (iRow in seq_len(nrow(normal_db))) {
      wpm[iRow, "WPM_Rank"] <- prod(sapply(normal_db[iRow, 
      ], as.numeric))
    }
    wpm_db <- wpm[, c("Alternative", "WPM_Rank")]
    return(wpm_db)
  }, error = function(cond) {
    stop(paste("E[P]", cond))
  }, warning = function(cond) {
    if (grepl("NAs intro", cond)) {
      return("W[P] Error: Some non numeric - alike value was found")
    }
  })
}