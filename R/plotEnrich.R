#' Plot the results of `enrichSEMs`
#' 
#' Generates a cladogram clustered by SNP Effect Matrix similarity and 
#' a heatmap representing the adjusted p-value of the enrichment.
#'
#' @param e The resulting data.table from `enrichSEMs`
#' @param threshold The adjusted p-value threshold for coloring SEMs
#' @param semList A `SNPEffectMatrixCollection` object
#'
#' @return a `ggtree` object
#'
#' @export
plotEnrich <- \(e, semList, threshold = 0.05) {
  sm <- NULL
  ppms <- convertSEMsToPPMs(sems(semList))
  motifs <- lapply(seq_along(ppms), 
                   function(i) universalmotif::create_motif(ppms[[i]], 
                                                            name = names(ppms[i])))
  
  for (i in 1:length(motifs)) {
    motifs[[i]]["altname"] <- semData(semList)[motifs[[i]]["name"], 
                                               "transcription_factor"] |>
      unlist()
  }
  
  comparisons <- universalmotif::compare_motifs(motifs, 
                                                method = "WPCC", 
                                                min.mean.ic = 0)
  labels <- lapply(motifs, function(x) x["altname"]) |> unlist()
  colnames(comparisons) <- labels
  rownames(comparisons) <- labels
  
  comparisons <- 1 - comparisons
  comparisons <- stats::as.dist(comparisons)
  
  comparisons <- ape::as.phylo(stats::hclust(comparisons))
  
  # significant motifs
  sig_motifs <- rownames(e)[e$padj <= 0.01]
  
  comparisons <- ggtree::groupOTU(comparisons,
                                  unlist(semData(semList)[sig_motifs, 
                                                          "transcription_factor"]),
                                  group_name = "sm")
  
  tree <- ggtree::ggtree(comparisons, layout = "fan", open.angle=35) +
    theme(legend.position = "", legend.title = element_blank()) +
    theme(plot.margin = margin(10, 10, 10, 100)) + 
    ggtree::geom_tiplab(aes(colour = sm), as_ylab = F, size = 1.5, offset = .05) +
    scale_color_manual(values=c("grey", "purple4"))
  
  heatmapData <- e[, "padj"] |> as.data.frame()
  colnames(heatmapData) <- "-log10(pval)"
  rownames(heatmapData) <- labels
  
  plt <- ggtree::gheatmap(tree, -log10(heatmapData), offset=-0.014, width=0.12,
                          colnames_angle=95, colnames_offset_y = 0.5, font.size = 3)
  
  plt <- plt +
    scale_fill_viridis_c(na.value = "white", direction = -1)
  print(plt)
  return(plt)
}
