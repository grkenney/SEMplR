#' Plot SEM scores for each nucleic acid in each position of the motif
#'
#' @param e data.table produced by running `SEMplR::enrichment()`
#'
#' @import ggplot2
#'
#' @return a `ggplot` with odds ratio and confidence intervals plotted for each
#' motif
#'
#' @export
#' 
plotEnrichVolcano <- function(e) {
  enrich_vol <- ggplot(e, aes(x = odds.ratio, 
                               y = -log10(adj.pvalue))) +
    geom_point(alpha=0.1) +
    geom_jitter(alpha=0.1) +
    theme_classic() +
    ylab("-Log10(adj.pval)") + xlab("Odds Ratio")
  
  return(enrich_vol)
}
