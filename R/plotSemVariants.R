#' Plot non-risk versus risk binding propensity for a single motif
#'
#' @param sempl_obj a SemplR object with scores populated
#' @param sem_id numeric, index of the variant within the SEMplR object to plot
#' @param label column in scores slot of sempl_obj to use for point labels
#'
#' @import ggplot2
#'
#' @return a `ggplot` scatter plot of non-risk binding propensity versus risk
#' binding propensity
#'
#' @export
plotSemVariants <- function(sempl_obj, sem_id, label = "varId") {
  nonRiskNorm <- riskNorm <- NA
  
  dt <- motifScores(sempl_obj, sem_id)
  
  var_plot <- ggplot2::ggplot(data = dt,
                              aes(x = nonRiskNorm, y = riskNorm)) +
    geom_vline(xintercept = 0, linetype = 2, col = 'grey') +
    geom_hline(yintercept = 0, linetype = 2, col  = 'grey') +
    geom_point(data = subset(dt, riskNorm < 0 & nonRiskNorm < 0),
               size = 1,
               color = 'grey') +
    geom_point(data = subset(dt, riskNorm > 0 & nonRiskNorm < 0),
               size = 1,
               color = 'firebrick') +
    geom_point(data = subset(dt, riskNorm < 0 & nonRiskNorm > 0),
               size = 1,
               color = '#1D91C0') +
    ggrepel::geom_text_repel(data = subset(dt, riskNorm > 0 & nonRiskNorm < 0),
                             mapping = aes(label = .data[[label]]),
                             size = 4,
                             color = 'firebrick') +
    ggrepel::geom_text_repel(data = subset(dt, riskNorm < 0 & nonRiskNorm > 0),
                             mapping = aes(label = .data[[label]]),
                             size = 4,
                             color = '#1D91C0') +
    geom_text(data = subset(dt, riskNorm > 0 & nonRiskNorm > 0),
              mapping = aes(label = .data[[label]]),
              size = 4,
              color = '#b05bc5') +
    scale_x_continuous(breaks = scales::pretty_breaks(),
                       limits = \(x) ifelse(abs(x) < 1, c(-1,1), x)) +
    scale_y_continuous(breaks = scales::pretty_breaks(),
                       limits = \(x) ifelse(abs(x) < 1, c(-1,1), x)) +
    labs(x = "Non-risk binding propensity",
         y = "Risk binding propensity") +
    theme_classic() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 18),
          axis.text = element_text(size=15))
  
  return(var_plot)
}
