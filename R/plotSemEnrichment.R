# Add columns with labels for plotting for significant transcription factors
prepPlotSemEnrichmentDf <- function(e, lab, sigThreshold, orThreshold) {
  e$lab <- lapply(1:nrow(e), 
                  function(i) 
                    ifelse((e$adj.pvalue[i] <= sigThreshold) &
                             (e$odds.ratio[i] > orThreshold),
                           as.character(e[[lab]][i]), ""))
  
  e$sig <- ifelse((e$adj.pvalue <= sigThreshold) & 
                    (e$odds.ratio > orThreshold), 
                  paste0("adj. p-value \u2264 ", sigThreshold), 
                  paste0("> ", sigThreshold))
  e$sig <- factor(e$sig, levels = c(paste0("adj. p-value \u2264 ", sigThreshold), 
                                    paste0("> ", sigThreshold)))
  return(e)
}

#' Plot motif enrichment results
#'
#' @param e SEMplR enrichment results data.table
#' @param lab column in enrichment results to label significant points by. 
#' Default is "semId".
#' @param sigThreshold adjusted pvalue threshold for coloring and labeling 
#' points. Default is 0.05
#' @param orThreshold odds ratio threshold for coloring and labeling points. 
#' Default is 1.25
#' @param sigCol color used to plot significant points
#'
#' @import ggplot2
#'
#' @return a `ggplot` with sem scores for each nucleic acid per position
#'
#' @export
plotSemEnrichment <- function(e, 
                              lab = "semId", 
                              sigThreshold = 0.05, 
                              orThreshold = 1.25,
                              sigCol = "dodgerblue2") {
  odds.ratio <- sig <- varId <- NULL
  
  e <- prepPlotSemEnrichmentDf(e, lab, sigThreshold, orThreshold)
  
  p <- ggplot(e, aes(x = stats::reorder(semId, odds.ratio), 
                     y = odds.ratio, col=sig)) + 
    geom_point() +
    theme_classic() +
    ggrepel::geom_text_repel(aes(semId, odds.ratio, label = lab), 
                             size = 3, show.legend = F, max.overlaps = Inf) +
    guides(colour  = guide_legend(position = "inside")) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size=10),
          legend.position = c(0.15, 0.95),
          legend.background = element_rect(fill=alpha('white', 0)),
          legend.key = element_rect(fill = "transparent", colour = "transparent")) +
    scale_color_manual(values = c(sigCol, "grey"),
                       breaks = c(paste0("adj. p-value \u2264 ", sigThreshold))) +
    labs(x = "Rank", y = "Odds Ratio", color='Adjusted\np-value') +
    coord_cartesian(xlim=c(0, nrow(e) + nrow(e)/8))
  
  return(p)
}
