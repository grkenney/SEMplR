#' Plot non-risk versus risk binding propensity
#'
#' @param sempl_obj a SemplR object with scores populated
#' @param variant_i numeric, index of the variant within the SEMplR object to plot
#'
#' @import ggplot2
#'
#' @return a `ggplot` scatter plot of non-risk binding propensity versus risk
#' binding propensity
#'
#' @export
plotSemMotifs <- \(sempl_obj, variant_i) {

  dt <- sempl_obj@scores[(variant_i*211-210):(variant_i*211), ]

  nonRiskNorm <- riskNorm <- sem <- NULL

  sem_motif_plot <- ggplot2::ggplot(data = dt,
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
                    mapping = aes(label = sem),
                    size = 4,
                    color = 'firebrick') +
    ggrepel::geom_text_repel(data = subset(dt, riskNorm < 0 & nonRiskNorm > 0),
                    mapping = aes(label = sem),
                    size = 4,
                    color = '#1D91C0') +
    geom_text(data = subset(dt, riskNorm > 0 & nonRiskNorm > 0),
              mapping = aes(label = sem),
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

  return(sem_motif_plot)
}
