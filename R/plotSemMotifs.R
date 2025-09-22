# validate label, variant, and cols parameters
.validatePlotSemMotifsInputs <- \(s, label, variant, cols) {
  semId <- varId <- NA
  # check that sem label is in the sem meta data
  if (!(label %in% colnames(semData(s)))) {
    rlang::abort("label not found in colnames(semData(s)).")
  }
  
  # check that variant is a valid id in s
  if (!(variant %in% scores(s)[, varId])) {
    rlang::abort(paste0("variant not found in SemplScores object. ",
                        variant, " is not in scores(s)[, varId]"))
  }
  
  # check that cols is length 2
  if (length(cols) != 2) {
    rlang::abort("cols must be a vector of length 2")
  }
}


#' Plot non-alt versus alt binding propensity for a single variant
#'
#' @param s a SemplScores object with scores populated
#' @param variant variant id to plot
#' @param label column in sem_metadata slot of semplObj to use for point labels
#' @param cols vector of length 2 with colors to use for plotting gained
#' and lost motifs respectively
#' 
#' @import ggplot2
#'
#' @return a `ggplot` scatter plot of non-alt binding propensity versus alt
#' binding propensity
#'
#' @export
#' 
#' @examples
#' library(VariantAnnotation)
#' data(sc)
#'
#' # create an SNP Effect Matrix (SEM)
#' sem <- matrix(rnorm(12), ncol = 4)
#' colnames(sem) <- c("A", "C", "G", "T")
#' 
#' # create a VRanges object
#' vr <- VRanges(seqnames = "chr12",
#'               ranges = 94136009, 
#'               ref = "G", alt = "C")
#' 
#' # calculate binding propensity
#' s <- scoreVariants(vr, sc, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
#' 
#' plotSemMotifs(s, "chr12:94136009:G>C", label = "transcription_factor")
#' 
plotSemMotifs <- \(s, variant, label = "semId",
                   cols = c("#F8766D", "dodgerblue2")) {
  refNorm <- altNorm <- varId <- sem <- NULL
  .validatePlotSemMotifsInputs(s = s, label = label, 
                               variant = variant, cols = cols)
  
  ix <- variant == scores(s)[, varId]
  dt <- scores(s)[ix, ]

  dt <- merge(dt, semData(s), 
              by.x = "semId", by.y = data.table::key(semData(s)))

  sem_motif_plot <- ggplot2::ggplot(data = dt,
                                    aes(x = refNorm, y = altNorm)) +
    geom_vline(xintercept = 0, linetype = 2, col = 'grey') +
    geom_hline(yintercept = 0, linetype = 2, col  = 'grey') +
    geom_point(data = subset(dt, altNorm < 0 & refNorm < 0),
               size = 1,
               color = 'grey') +
    geom_point(data = subset(dt, altNorm > 0 & refNorm < 0),
               size = 1,
               color = cols[1]) +
    geom_point(data = subset(dt, altNorm < 0 & refNorm > 0),
               size = 1,
               color = cols[2]) +
    geom_point(data = subset(dt, altNorm > 0 & refNorm > 0),
               size = 1,
               color = 'grey') +
    ggrepel::geom_text_repel(data = subset(dt, altNorm > 0 & refNorm < 0),
                    mapping = aes(label = .data[[label]]),
                    size = 4,
                    color = cols[1]) +
    ggrepel::geom_text_repel(data = subset(dt, altNorm < 0 & refNorm > 0),
                    mapping = aes(label = .data[[label]]),
                    size = 4,
                    color = cols[2]) +
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
