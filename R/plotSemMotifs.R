plotSemMotifs <- \(dt) {
  library(ggplot2)
  ggplot(data = dt, aes(x = nonrisk, y = risk)) +
    geom_vline(xintercept = 0, linetype = 2, col = 'grey') +
    geom_hline(yintercept = 0, linetype = 2, col  = 'grey') +
    geom_point(data = subset(dt, risk < 0 & nonrisk < 0),
               size = 1,
               color = 'grey') +
    geom_text(data = subset(dt, risk > 0 & nonrisk < 0),
              mapping = aes(label = motif),
              size = 4,
              color = 'firebrick') +
    geom_text(data = subset(dt, risk < 0 & nonrisk > 0),
              mapping = aes(label = motif),
              size = 4,
              color = '#1D91C0') +
    geom_text(data = subset(dt, risk > 0 & nonrisk > 0),
              mapping = aes(label = motif),
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
}
