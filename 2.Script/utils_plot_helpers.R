# ----------------------------------------------------------------------------
# Helper | generate_plot_data(): builds the GESECA pathway heatmap matrix
#
# Paper    : Ligand-induced ubiquitination unleashes LAG3 immune checkpoint
#            function by hindering membrane sequestration of signaling motifs
# Original : 2.Script/plotInfo.r
# Note     : Sourced by 09_functional_enrichment_geseca.R
# ----------------------------------------------------------------------------

generate_plot_data <- function(gesecaRes, pathways, E, center = TRUE, scale = FALSE, titles = colnames(E)) {
  # Filter gesecaRes and pathways
  gesecaRes <- gesecaRes[gesecaRes$pathway %in% names(pathways), ]
  pathways <- pathways[gesecaRes$pathway]
  
  # Scale the expression data
  E <- t(base::scale(t(E), center = center, scale = scale))
  colnames(E) <- titles
  
  # Process pathways
  pathways <- lapply(pathways, function(p) {
    unname(as.vector(na.omit(match(p, rownames(E)))))
  })
  
  # Calculate projections
  prjs <- t(do.call(cbind, lapply(pathways, function(p) {
    scale(colSums(E[p, , drop = FALSE]))
  })))
  rownames(prjs) <- names(pathways)
  
  # Convert to long format correctly
  prjspd <- data.frame(
    pathway = rep(rownames(prjs), each = ncol(prjs)),
    sample = rep(colnames(prjs), times = nrow(prjs)),
    value = as.vector(t(prjs))
  )
  
  # Add factor levels
  prjspd$pathway <- factor(prjspd$pathway, levels = rev(rownames(prjs)))
  prjspd$sample <- factor(prjspd$sample, levels = titles)
  
  # Add annotation data
  annotation <- gesecaRes[match(unique(prjspd$pathway), gesecaRes$pathway), ]
  
  list(
    plot_data = prjspd,
    annotation = annotation
  )
}