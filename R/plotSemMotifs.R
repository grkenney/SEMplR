#' Plot non-alt versus alt binding propensity for a single variant
#'
#' @param semplObj a SemplScores object with scores populated
#' @param variant variant id to plot
#' @param label column in sem_metadata slot of semplObj to use for point labels
#'
#' @import ggplot2
#'
#' @return a `ggplot` scatter plot of non-alt binding propensity versus alt
#' binding propensity
#'
#' @export
plotSemMotifs <- \(semplObj, variant, label = "semId") {
  refNorm <- altNorm <- varId <- sem <- NULL
  
  if (variant %in% scores(semplObj)[, varId]) {
    dt <- scores(semplObj)[varId == variant]
  } else {
    stop(paste0("variant not found in SemplScores object. ",
                variant, " is not in scores(semplObj)[, varId]"))
  }
  
  dt <- merge(dt, semData(semplObj), 
              by.x = "semId", by.y = data.table::key(semData(semplObj)))

  sem_motif_plot <- ggplot2::ggplot(data = dt,
                                    aes(x = refNorm, y = altNorm)) +
    geom_vline(xintercept = 0, linetype = 2, col = 'grey') +
    geom_hline(yintercept = 0, linetype = 2, col  = 'grey') +
    geom_point(data = subset(dt, altNorm < 0 & refNorm < 0),
               size = 1,
               color = 'grey') +
    geom_point(data = subset(dt, altNorm > 0 & refNorm < 0),
               size = 1,
               color = 'firebrick') +
    geom_point(data = subset(dt, altNorm < 0 & refNorm > 0),
               size = 1,
               color = '#1D91C0') +
    geom_point(data = subset(dt, altNorm > 0 & refNorm > 0),
               size = 1,
               color = '#b05bc5') +
    ggrepel::geom_text_repel(data = subset(dt, altNorm > 0 & refNorm < 0),
                    mapping = aes(label = .data[[label]]),
                    size = 4,
                    color = 'firebrick') +
    ggrepel::geom_text_repel(data = subset(dt, altNorm < 0 & refNorm > 0),
                    mapping = aes(label = .data[[label]]),
                    size = 4,
                    color = '#1D91C0') +
    ggrepel::geom_text_repel(data = subset(dt, altNorm > 0 & refNorm > 0),
                             mapping = aes(label = .data[[label]]),
                             size = 4,
                             color = '#b05bc5') +
    scale_x_continuous(breaks = scales::pretty_breaks(),
                       limits = \(x) ifelse(abs(x) < 1, c(-1,1), x)) +
    scale_y_continuous(breaks = scales::pretty_breaks(),
                       limits = \(x) ifelse(abs(x) < 1, c(-1,1), x)) +
    labs(x = "ref binding propensity",
         y = "alt binding propensity") +
    theme_classic() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 18),
          axis.text = element_text(size=15))

  return(sem_motif_plot)
}
