#' Plot SEM scores for each nucleic acid in each position of the motif
#'
#' @param e data.table produced by running `SEMplR::enrichment()`
#' @param orc (optional) odds ratio cut-off to filter enrichment results. Only 
#' results greater than the cut-off will be included.
#' @param pvc (optional) adjusted p-value cut-off to filter enrichment results.
#' Only results less than the cut-off will be included.
#'
#' @import ggplot2
#'
#' @return a `ggplot` with odds ratio and confidence intervals plotted for each
#' motif
#'
#' @export
#' 
plotOddsRatio <- function(e, orc = NULL, pvc = NULL) {
  if (!is.null(orc)) {
    e <- e[odds.ratio > orc]
  }
  
  if (!is.null(pvc)) {
    e <- e[adj.pvalue < pvc]
  }
  
  enrich_plot <- ggplot(e, aes(x = reorder(motif, odds.ratio), 
                               y = odds.ratio)) +
    geom_point() +
    geom_errorbar(aes(ymin=ci.lower, ymax=ci.upper)) +
    theme_classic() +
    coord_flip() +
    xlab("SEM") + ylab("Odds Ratio")
  
  return(enrich_plot)
}
