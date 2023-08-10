compute_rotated_pca <- function(.data = NULL, naxis = 2) {
  .data <- na.omit(.data)
  pca <- ade4::dudi.pca(as.data.frame(.data),
    scannf = FALSE, nf = naxis, center = TRUE, scale = TRUE
  )
  pca_rotated <- psych::principal(.data,
    rotate = "varimax", nfactors = naxis,
    scores = TRUE
  )

  return(list(normal = pca, rotated = pca_rotated))
}

plot_rotated_pca <- function(pca_rotated = NULL, axis = c(1, 2)) {
  {
    adegraphics::adegpar()$plabels
    ade4::s.corcircle(
      pca_rotated$rotated$loadings[1:nrow(pca_rotated$rotated$loadings), ],
      xax = axis[1],
      yax = axis[2], box = TRUE
    )
    mtext(paste0("RC ", axis[1]), side = 1, adj = .5, line = -1)
    mtext(paste0("RC ", axis[2]), side = 2, adj = .5, line = -5)
    rotated_plot <- grDevices::recordPlot()
  }
  pca_plot <- factoextra::fviz_pca_var(pca_rotated$normal,
    axes = axis,
    col.var = "contrib", # Color by contributions to the PC
    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
    repel = TRUE # Avoid text overlapping
  )
  return(list(rotated = rotated_plot, normal = pca_plot))
}
