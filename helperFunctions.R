parse_flowjo <- function(files,
                         wsp_file,
                         group = "All Samples",
                         plot = FALSE) {
  wsp <- flowWorkspace::openWorkspace(wsp_file)
  o <- capture.output(
    gates <- suppressMessages(flowWorkspace::parseWorkspace(wsp, group))
  )
  files_in_wsp <- gates@data@origSampleVector
  counts <- as.numeric(gsub(".*_([0-9]*)$", "\\1", files_in_wsp))
  result <- list()
  for(file in files){
    print(paste0("Processing ", file))
    file_id <- grep(gsub(".*/", "", file), files_in_wsp)
    if(length(file_id) == 0) {stop("File not found. Files available: ",
                                   gsub("_[0-9]*$", "\n", files_in_wsp))}
    gate_names <- flowWorkspace::getNodes(gates[[file_id]], path = "auto")
    gatingMatrix <- matrix(FALSE,
                           nrow = counts[file_id],
                           ncol = length(gate_names),
                           dimnames = list(NULL, gate_names))
    for (gate in gate_names) {
      gatingMatrix[, gate] <- flowWorkspace::getIndiceMat(gates[[file_id]],
                                                          gate)
    }
    ff <- flowWorkspace::getData(gates[[file_id]], "root")
    ff@exprs[, "Time"] <- ff@exprs[, "Time"] * 100
    result[[file]] <- list("flowFrame" = ff,
                           "gates" = gatingMatrix)
    
    if (plot) {
      flowWorkspace::plot(gates[[file_id]])
    }
  }
  if (length(files) == 1){
    result <- result[[1]]
  } else {
    result <- list(flowSet = flowCore::flowSet(lapply(result, function(x) x$flowFrame)),
                   gates = lapply(result, function(x) x$gates))
  }
  return(result)
}

manual_vector <- function(manual_matrix, cell_types){
  
  if(is.list(manual_matrix)){ manual_matrix <- do.call(rbind, manual_matrix) }
  
  manual <- rep("Unknown",nrow(manual_matrix))
  for(cellType in cell_types){
    manual[manual_matrix[,cellType]] <- cellType
  }
  manual <- factor(manual, levels=c("Unknown",cell_types))
  return(manual)
}

label_metaclusters <- function(fsom, manual_labels){
  counts <- as.matrix(table(FlowSOM::GetMetaclusters(fsom),
                            manual_labels))
  metacluster_names <- apply(counts,
                             1,
                             function(x) colnames(counts)[which.max(x)])
  metacluster_names <- number_duplicates(metacluster_names)
  return(metacluster_names)
}

label_clusters <- function(fsom, manual_labels){
  counts <- as.matrix(table(FlowSOM::GetClusters(fsom),
                            manual_labels))
  cluster_names <- apply(counts,
                         1,
                         function(x) colnames(counts)[which.max(x)])
  return(cluster_names)
}

number_duplicates <- function(x){
  counts <- table(x)
  for (value in names(counts)) {
    if (counts[value] > 1) {
      x[which(x == value)] <- paste0(value, "_", seq_len(counts[value]))
    }
  }
  return(x)
}

PlotDiscreteVariable <- function (fsom, 
                                  variable, 
                                  view = "MST", 
                                  main = NULL, 
                                  colors = grDevices::colorRampPalette(c(
                                    "#00007F", "blue", "#007FFF", "cyan", 
                                    "#7FFF7F", "yellow", "#FF7F00", "red", 
                                    "#7F0000"))(length(levels(variable))), 
                                  symmetric = FALSE, 
                                  lim = NULL, 
                                  backgroundValues = NULL, 
                                  backgroundColor = function(n) {
                                    grDevices::rainbow(n, alpha = 0.3)
                                  }, 
                                  backgroundLim = NULL, 
                                  backgroundBreaks = NULL) {
  switch(view, 
         MST = {
           layout <- fsom$MST$l
           lty <- 1
         }, grid = {
           layout <- as.matrix(fsom$map$grid)
           lty <- 0
         }, tSNE = {
           layout <- fsom$MST$l2
           lty <- 0
         }, stop("The view should be MST, grid or tSNE. tSNE will only work\n                   
          if you specified this when building the MST."))
  
  if (!is.null(backgroundValues)) {
    background <- computeBackgroundColor(backgroundValues, 
                                         backgroundColor, backgroundLim, backgroundBreaks)
  }
  else {
    background <- NULL
  }
  
  oldpar <- graphics::par(no.readonly = TRUE)
  graphics::par(mar = c(1, 1, 1, 1))
  graphics::layout(matrix(c(1, 2), 1, 2, byrow = TRUE), widths = c(1, 
                                                                   2), heights = c(1))
  graphics::plot.new()
  legend("center", legend = levels(variable), fill = colors[levels(variable)], 
         cex = 0.7, ncol = 2, bty = "n")
  
  f <- fsom
  igraph::V(f$MST$graph)$color <- colors[as.character(variable)]
  #colorPalette(100)[as.numeric(cut(c(lim, 
  #                                   variable), breaks = 100))[-c(1, 2)]]
  igraph::plot.igraph(f$MST$graph, 
                      layout = layout, 
                      vertex.size = fsom$MST$size, 
                      vertex.label = NA, 
                      main = main, 
                      edge.lty = lty, 
                      mark.groups = background$groups, 
                      mark.col = background$col[background$values], 
                      mark.border = background$col[background$values])
  graphics::layout(1)
  graphics::par(oldpar)
}

# Variant of the function in the FlowSOM package to adapt the size of the
# background circles.
PlotPies <- PlotPies <- function (fsom, 
                                  cellTypes, 
                                  view = "MST", 
                                  colorPalette = grDevices::colorRampPalette(
                                    c("white", "#00007F", "blue", "#007FFF", 
                                      "cyan", "#7FFF7F", "yellow", "#FF7F00", 
                                      "red")), 
                                  backgroundValues = NULL, 
                                  backgroundColor = function(n) {
                                    grDevices::rainbow(n, alpha = 0.3)}, 
                                  backgroundLim = NULL, 
                                  backgroundBreaks = NULL, 
                                  legend = TRUE, 
                                  main = "") 
{
  if (!is.factor(cellTypes)) {
    cellTypes <- as.factor(cellTypes)
  }
  t <- table(factor(fsom$map$mapping[, 1], levels = seq_along(fsom$MST$size)), 
             factor(cellTypes, levels = c(levels(cellTypes), "empty")))
  t[rowSums(t) == 0, "empty"] <- 1
  data <- unlist(apply(t, 1, list), recursive = FALSE)
  colors <- list(c(colorPalette(length(levels(cellTypes))), 
                   "#000000"))
  colors <- rep(colors, length(fsom$MST$size))
  switch(view, MST = {
    layout <- fsom$MST$l
    lty <- 1
  }, grid = {
    layout <- as.matrix(fsom$map$grid)
    lty <- 0
  }, tSNE = {
    layout <- fsom$MST$l2
    lty <- 0
  }, stop("The view should be MST, grid or tSNE. tSNE will only work\n                   
          if you specified this when building the MST."))
  if (!is.null(backgroundValues)) {
    background <- FlowSOM:::computeBackgroundColor(backgroundValues, 
                                                   backgroundColor, backgroundLim, backgroundBreaks)
    backgroundSize <- fsom$MST$size
    backgroundSize[backgroundSize == 0] <- 3
  }
  else {
    background <- NULL
    backgroundSize <- NULL
  }
  oldpar <- graphics::par(no.readonly = TRUE)
  graphics::par(mar = c(1, 1, 1, 1))
  if (legend) {
    if (!is.null(backgroundValues)) {
      graphics::layout(matrix(c(1, 3, 2, 3), 2, 2, byrow = TRUE), 
                       widths = c(1, 2), heights = c(1))
    }
    else {
      graphics::layout(matrix(c(1, 2), 1, 2, byrow = TRUE), 
                       widths = c(1, 2), heights = c(1))
    }
    graphics::plot.new()
    legend("center", legend = levels(cellTypes), fill = colors[[1]], 
           cex = 0.7, ncol = 1, bty = "n")
    if (!is.null(backgroundValues)) {
      FlowSOM:::PlotBackgroundLegend(backgroundValues, background)
    }
  }
  igraph::plot.igraph(fsom$MST$g, vertex.shape = "pie", vertex.label = NA, 
                      vertex.size = fsom$MST$size, vertex.pie = data, 
                      vertex.pie.color = colors, 
                      layout = layout, edge.lty = lty, 
                      mark.groups = background$groups, 
                      mark.col = background$col[background$values], 
                      mark.border = background$col[background$values], 
                      mark.expand = backgroundSize,
                      main = main)
  graphics::par(oldpar)
  graphics::layout(1)
}
