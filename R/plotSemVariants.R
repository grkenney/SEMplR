#' Plot non-alt versus alt binding propensity for a single motif
#'
#' @param semplObj a SemplScores object with scores populated
#' @param semId a single character vector matching a semId in the semplObj
#' @param label column in scores slot of semplObj to use for point labels
#' @param changedCols vector of length 2 with colors to use for plotting gained
#' and lost motifs respectively
#'
#' @import ggplot2
#'
#' @return a `ggplot` scatter plot of non-alt binding propensity versus alt
#' binding propensity
#'
#' @export
plotSemVariants <- function(semplObj, semId, label = "varId",
                            changedCols = c("#F8766D", "dodgerblue2")) {
  refNorm <- altNorm <- NA
  
  if (semId %in% scores(semplObj)[, semId]) {
    ix <- semId == scores(semplObj)[, semId]
    dt <- scores(semplObj)[ix, ]
  } else {
    stop(paste0("semId not found in SemplScores object. ",
                semId, " is not in scores(semplObj)[, semId]"))
  }
  
  if (length(changedCols) != 2) {
    stop("changedCols must be a vector of length 2")
  }
  
  var_plot <- ggplot2::ggplot(data = dt,
                              aes(x = refNorm, y = altNorm)) +
    geom_vline(xintercept = 0, linetype = 2, col = 'grey') +
    geom_hline(yintercept = 0, linetype = 2, col  = 'grey') +
    geom_point(data = subset(dt, altNorm < 0 & refNorm < 0),
               size = 1,
               color = 'grey') +
    geom_point(data = subset(dt, altNorm > 0 & refNorm < 0),
               size = 1,
               color = changedCols[1]) +
    geom_point(data = subset(dt, altNorm < 0 & refNorm > 0),
               size = 1,
               color = changedCols[2]) +
    geom_point(data = subset(dt, altNorm > 0 & refNorm > 0),
               size = 1,
               color = 'grey') +
    ggrepel::geom_text_repel(data = subset(dt, altNorm > 0 & refNorm < 0),
                             mapping = aes(label = .data[[label]]),
                             size = 4,
                             color = changedCols[1]) +
    ggrepel::geom_text_repel(data = subset(dt, altNorm < 0 & refNorm > 0),
                             mapping = aes(label = .data[[label]]),
                             size = 4,
                             color = changedCols[2]) +
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
  
  return(var_plot)
}
