DimHeatmap_custom<-function (object, dims = 1, nfeatures = 30, cells = NULL, reduction = "pca", 
          disp.min = -2.5, disp.max = NULL, balanced = TRUE, projected = FALSE, 
          ncol = NULL, fast = TRUE, raster = TRUE, slot = "scale.data", 
          assays = NULL, combine = TRUE, s.genes, g2m.genes) 
{
  ncol <- ncol %||% ifelse(test = length(x = dims) > 2, yes = 3, 
                           no = length(x = dims))
  plots <- vector(mode = "list", length = length(x = dims))
  assays <- assays %||% DefaultAssay(object = object)
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                   yes = 2.5, no = 6)
  if (!DefaultAssay(object = object[[reduction]]) %in% assays) {
    warning("The original assay that the reduction was computed on is different than the assay specified")
  }
  cells <- cells %||% ncol(x = object)
  if (is.numeric(x = cells)) {
    cells <- lapply(X = dims, FUN = function(x) {
      cells <- TopCells(object = object[[reduction]], dim = x, 
                        ncells = cells, balanced = balanced)
      if (balanced) {
        cells$negative <- rev(x = cells$negative)
      }
      cells <- unlist(x = unname(obj = cells))
      return(cells)
    })
  }
  if (!is.list(x = cells)) {
    cells <- lapply(X = 1:length(x = dims), FUN = function(x) {
      return(cells)
    })
  }
  features <- lapply(X = dims, FUN = TopFeatures, object = object[[reduction]], 
                     nfeatures = nfeatures, balanced = balanced, projected = projected)
  features.all <- unique(x = unlist(x = features))
  if (length(x = assays) > 1) {
    features.keyed <- lapply(X = assays, FUN = function(assay) {
      features <- features.all[features.all %in% rownames(x = object[[assay]])]
      if (length(x = features) > 0) {
        return(paste0(Key(object = object[[assay]]), 
                      features))
      }
    })
    features.keyed <- Filter(f = Negate(f = is.null), x = features.keyed)
    features.keyed <- unlist(x = features.keyed)
  }else {
    features.keyed <- features.all
    DefaultAssay(object = object) <- assays
  }
  data.all <- FetchData(object = object, vars = features.keyed, 
                        cells = unique(x = unlist(x = cells)), slot = slot)
  data.all <- MinMax(data = data.all, min = disp.min, max = disp.max)
  data.limits <- c(min(data.all), max(data.all))
  if (fast) {
    nrow <- floor(x = length(x = dims)/3.01) + 1
    orig.par <- par()$mfrow
    par(mfrow = c(nrow, ncol))
  }
  for (i in 1:length(x = dims)) {
    dim.features <- c(features[[i]][[2]], rev(x = features[[i]][[1]]))
    dim.features <- rev(x = unlist(x = lapply(X = dim.features, 
                                              FUN = function(feat) {
                                                return(grep(pattern = paste0(feat, "$"), x = features.keyed, 
                                                            value = TRUE))
                                              })))
    dim.cells <- cells[[i]]
    data.plot <- data.all[dim.cells, dim.features]
    
    names(data.plot) = ifelse(names(data.plot)%in%s.genes, paste0(names(data.plot), '(S)'), names(data.plot))
    names(data.plot) = ifelse(names(data.plot)%in%g2m.genes, paste0(names(data.plot), '(G2M)'), names(data.plot))
    
    if (fast) {
      SingleImageMap(data = data.plot, title = paste0(Key(object = object[[reduction]]), 
                                                      dims[i]), order = dim.cells)
    }
    else {
      plots[[i]] <- SingleRasterMap(data = data.plot, raster = raster, 
                                    limits = data.limits, cell.order = dim.cells, 
                                    feature.order = dim.features)
    }
  }
  if (fast) {
    par(mfrow = orig.par)
    return(invisible(x = NULL))
  }
  if (combine) {
    plots <- wrap_plots(plots, ncol = ncol, guides = "collect")
  }
  return(plots)
}
