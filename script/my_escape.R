

.pull.Enrich <- function(sc, enrichment.name) {
  if (inherits(sc, "Seurat")) {
    values <- t(sc[[enrichment.name]]["data"])
  } else if (inherits(sc, "SingleCellExperiment")) {
    if(length(assays(altExp(sc))) == 1) {
      values <- t(assay(altExps(sc)[[enrichment.name]]))
    }
  }
}

.GS.check <- function(gene.sets) {
  if(is.null(gene.sets)) {
    stop("Please provide the gene.sets you would like to use for 
            the enrichment analysis")
  }
  egc <- gene.sets
  if(inherits(egc, what = "GeneSetCollection")){
    egc <- GSEABase::geneIds(egc) # will return a simple list, 
    #which will work if a matrix is supplied to GSVA
  }
  return(egc)
}

.cntEval <- function(obj, 
                     assay = "RNA", 
                     type = "counts") {
  if (inherits(x = obj, what = "Seurat")) {
    cnts <- obj@assays[[assay]][type]
  } else if (inherits(x = obj, what = "SingleCellExperiment")) {
    pos <- ifelse(assay == "RNA", "counts", assay) 
    if(assay == "RNA") {
      cnts <- assay(obj,pos)
    } else {
      cnts <- assay(altExp(obj), pos)
    }
  } else {
    cnts <- obj
  }
  cnts <- cnts[rowSums2(cnts) != 0,]
  return(cnts)
}

.cntEval_mod <- function(obj, 
                         assay = "RNA", 
                         type = "counts") {
  print("check_seurat")
  if (inherits(x = obj, what = "Seurat")) {
    cnts <- obj@assays[[assay]][type]
  } else if (inherits(x = obj, what = "SingleCellExperiment")) {
    pos <- ifelse(assay == "RNA", "counts", assay) 
    if(assay == "RNA") {
      cnts <- assay(obj,pos)
    } else {
      cnts <- assay(altExp(obj), pos)
    }
  } else {
    cnts <- obj
  }
  print("calculate row_sum")
  cnts <- cnts[rowSums(cnts,na.rm = TRUE) != 0,] #modified this line by adding na.rm
  print("finish calculate row_sum")
  return(cnts)
}

is_seurat_or_se_object <- function(obj) {
  is_seurat_object(obj) || is_se_object(obj)
}
is_seurat_object <- function(obj) inherits(obj, "Seurat")
is_se_object <- function(obj) inherits(obj, "SummarizedExperiment")

.adding.Enrich <- function(sc, enrichment, enrichment.name) {
  if (inherits(sc, "Seurat")) {
    new.assay <- suppressWarnings(CreateAssayObject(
      data = as.matrix(t(enrichment))))
    
    sc[[enrichment.name]] <- new.assay
  } else if (inherits(sc, "SingleCellExperiment")) {
    altExp(sc, enrichment.name) <- SummarizedExperiment(assays = t(enrichment))
    names(assays(altExp(sc, enrichment.name))) <- enrichment.name
  }
  return(sc)
}










performNormalization_edited <- function (input.data, assay = NULL, gene.sets = NULL, make.positive = FALSE, 
                                         scale.factor = NULL) 
{
  
  library(pbapply)
  pbo = pboptions(type="txt") # !!!
  library(dplyr)
  if (is_seurat_or_se_object(input.data)) {
    print(".pull.Enrich")
    enriched <- .pull.Enrich(input.data,enrichment.name =  assay)
  }
  else {
    enriched <- input.data
  }
  if (!is.null(scale.factor) & length(scale.factor) != dim(input.data)[2]) {
    stop("If using a vector as a scale factor, please ensure the length matches the number of cells.")
  }
  print(".GS.check")
  egc <- .GS.check(gene.sets)
  names(egc) <- stringr::str_replace_all(names(egc), "_", "-")
  egc <- egc[names(egc) %in% colnames(enriched)]
  print(".cntEval_mod")
  cnts <- .cntEval_mod(input.data, assay = "RNA", type = "counts")
  if (is.null(scale.factor)) {
    print("Calculating features per cell... (1/3)")
    non_zero_indices <- pblapply(seq_len(ncol(cnts)), function(y) {
      which(cnts[, y] != 0)
    })
    print("Calculating features per cell... (2/3)")
    egc_indices <- pblapply(egc, function(x) {
      which(rownames(cnts) %in% x)
    })
    print("Calculating features per cell...(3/3)")
    egc.size <- pblapply(egc_indices, function(gene_set_indices) {
      sapply(non_zero_indices, function(sample_indices) {
        length(intersect(sample_indices, gene_set_indices))
      })
    })
  }

  print("Normalizing enrichment scores per cell...")
  normalized.values <- pblapply(seq_len(ncol(enriched)), function(x) {
    
    if (!is.null(scale.factor)) {
      enriched[, x] <- enriched[, x]/scale.factor
    }
    else {
      gene.set <- unlist(egc.size[colnames(enriched)[x]])
      if (any(gene.set == 0)) {
        gene.set[which(gene.set == 0)] <- 1
      }
      enriched[, x] <- enriched[, x]/gene.set
    }
    if (any(enriched[, x] < 0) & make.positive) {
      enriched[, x] <- enriched[, x] + abs(min(enriched[, 
                                                        x]))
    }
    enriched[, x]
  })
  
  
  # input.data = gbm.list_positive_TRUE; assay = "escape.ssGSEA"; gene.sets = gs_all; make.positive = TRUE; scale.factor = NULL
  #   for (x in seq_len(ncol(enriched))) {
  #     print(x)
  #   if (!is.null(scale.factor)) {
  #     enriched[, x] <- enriched[, x]/scale.factor
  #   }
  #   else {
  #     gene.set <- unlist(egc.size[colnames(enriched)[x]])
  #     if (any(gene.set == 0)) {
  #       gene.set[which(gene.set == 0)] <- 1
  #     }
  #     enriched[, x] <- enriched[, x]/gene.set
  #   }
  #   if (any(enriched[, x] < 0) & make.positive) {
  #     enriched[, x] <- enriched[, x] + abs(min(enriched[, 
  #                                                       x]))
  #   }
  # }
  
  
  normalized.enriched <- do.call(cbind, normalized.values)
  colnames(normalized.enriched) <- colnames(enriched)
  if (is_seurat_or_se_object(input.data)) {
    input.data <- .adding.Enrich(input.data, normalized.enriched, 
                                 paste0(assay, "_normalized"))
    return(input.data)
  }
  else {
    return(normalized.enriched)
  }
}



.prepData <- function(input.data, assay, gene.set, group.by, split.by, facet.by) {
  
  
  if (inherits(x=input.data, what ="Seurat") || 
      inherits(x=input.data, what ="SummarizedExperiment")) {
    enriched <- .makeDFfromSCO(input.data, assay, gene.set, group.by, split.by, facet.by)
    if(length(gene.set) == 1 && gene.set == "all") {
      gene.set <- colnames(enriched)[colnames(enriched) %!in% c(group.by, split.by, facet.by)]
      gene.set <- gene.set[!grepl("meta", gene.set)]
    }
  } else if (!is_seurat_or_se_object(input.data)) {
    if(length(gene.set) == 1 && gene.set == "all") {
      gene.set <- colnames(input.data)
      gene.set <- gene.set[gene.set %!in% c(group.by, split.by, facet.by)]
    } 
    enriched <- data.frame(input.data[,c(gene.set,group.by, split.by, facet.by)])
  }
  
  colnames(enriched) <- c(gene.set, group.by, split.by, facet.by)
  return(enriched)
}

.makeDFfromSCO <- function(input.data, 
                           assay = "escape", 
                           gene.set = NULL,
                           group.by = NULL, 
                           split.by = NULL, 
                           facet.by = NULL) {
  if(is.null(assay)){
    stop("Please add the assay name in which to plot from")
  }
  columns <- unique(c(group.by, split.by, facet.by))
  
  cnts <- .cntEval_mod(input.data,  #change this code to use modified version
                       assay = assay, 
                       type = "data")
  if(length(gene.set) == 1 && gene.set == "all") {
    gene.set <- rownames(cnts)
  }
  meta <- .grabMeta(input.data)
  
  
  ggVennDiagram::ggVennDiagram(list(
    gene.set = gene.set,
    rowname_cnts = rownames(cnts)
  ))
  
  if(length(gene.set) == 1) {
    enriched <- data.frame(cnts[gene.set,], meta[,columns])
  } else {
    enriched <- data.frame(t(cnts[gene.set,]), meta[,columns])
  }
  colnames(enriched) <- c(gene.set, columns)
  return(enriched)
}

.grabMeta <- function(sc) {
  if (is_seurat_object(sc)) {
    meta <- data.frame(sc[[]], slot(sc, "active.ident"))
    colnames(meta)[length(meta)] <- "ident"
    
  } else if (is_se_object(sc)){
    meta <- data.frame(colData(sc))
    rownames(meta) <- sc@colData@rownames
    clu <- which(colnames(meta) == "ident")
    colnames(meta)[clu] <- "ident"
  } else {
    stop("Object indicated is not of class 'Seurat' or 
            'SummarizedExperiment', make sure you are using
            the correct data.")
  }
  return(meta)
}

"%!in%" <- Negate("%in%")


heatmapEnrichment_modified <- function (input.data, assay = NULL, group.by = NULL, gene.set.use = "all", 
                                        cluster.rows = FALSE, cluster.columns = FALSE, scale = FALSE, 
                                        facet.by = NULL, palette = "inferno") 
{
  library(reshape2)
  library(ggplot2)
  options(dplyr.summarise.inform = FALSE)
  if (is.null(group.by)) {
    group.by <- "ident"
  }
  enriched <- .prepData(input.data, assay, gene.set.use, group.by, 
                        NULL, facet.by)
  if (length(gene.set.use) == 1 && gene.set.use == "all") {
    gene.set <- colnames(enriched)[colnames(enriched) %!in% 
                                     c(group.by, facet.by)]
  }
  else {
    gene.set <- gene.set.use
  }
  if (!is.null(facet.by)) {
    enriched.summary <- enriched %>% group_by(.data[[group.by]], 
                                              .data[[facet.by]]) %>% summarise(across(which(colnames(enriched) %in% 
                                                                                              gene.set), mean)) %>% as.data.frame()
  }
  else {
    enriched.summary <- enriched %>% group_by(.data[[group.by]]) %>% 
      summarise(across(which(colnames(enriched) %in% gene.set), 
                       mean)) %>% as.data.frame()
  }
  if (scale) {
    enriched.summary[, gene.set] <- apply(enriched.summary[, 
                                                           gene.set], 2, scale)
  }
  reformated.enriched <- suppressMessages(melt(enriched.summary))
  if (cluster.rows) {
    row.order <- gene.set[hclust(dist(t(enriched.summary[, 
                                                         gene.set]), method = "euclidean"), method = "ward.D2")$order]
    reformated.enriched[, "variable"] <- factor(reformated.enriched[, 
                                                                    "variable"], levels = row.order)
  }
  if (cluster.columns) {
    column.order <- unique(enriched.summary[, group.by][hclust(dist(enriched.summary[, 
                                                                                     gene.set], method = "euclidean"), method = "ward.D2")$order])
    reformated.enriched[, group.by] <- factor(reformated.enriched[, 
                                                                  group.by], levels = as.vector(column.order))
  }
  plot <- ggplot(reformated.enriched, mapping = aes(x = reformated.enriched[, 
                                                                            group.by], y = variable, fill = value)) + geom_tile(color = "black", 
                                                                                                                                linewidth = 0.5) + scale_y_discrete(expand = c(0, 0)) + 
    scale_x_discrete(expand = c(0, 0)) + labs(fill = "Enrichment Score") + 
    guides(fill = guide_colorbar(title.position = "top", 
                                 title.hjust = 0.5)) + coord_equal() + scale_fill_gradientn(colors = .colorizer(palette, 
                                                                                                                11)) + theme_classic() + theme(axis.title = element_blank(), 
                                                                                                                                               axis.ticks = element_blank(), legend.direction = "horizontal", 
                                                                                                                                               legend.position = "bottom",
                                                                                                                                               axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
                                                                                                                )
  if (!is.null(facet.by)) {
    plot <- plot + facet_grid(as.formula(paste(". ~", facet.by)))
  }
  return(plot)
}

.colorizer <- function(palette = "inferno", 
                       n= NULL) {
  colors <- hcl.colors(n=n, palette = palette, fixup = TRUE)
  return(colors)
}
