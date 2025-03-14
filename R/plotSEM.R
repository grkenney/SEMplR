#' Plot SEM scores for each nucleic acid in each position of the motif
#'
#' @param s A SNPEffectMatrix or SNPEffectMatrixCollection object.
#' @param motif sem id to plot. If providing a SNPEffectMatrixCollection object,
#' must match a name in the list. If providing a single SNPEffectMatrix object,
#' this parameter is ignored and the SNPEffectMatrix's semId is used for 
#' plotting.
#' @param motifSeq (optional) character sequence to color on plot
#' @param highlight (optional) index of variant location on reference
#' @param hcol (optional) color of motifSeq basepairs
#'
#' @import ggplot2
#'
#' @return a `ggplot` with sem scores for each nucleic acid per position
#'
#' @export
#' 
#' @examples
#' library(VariantAnnotation)
#'
#' # create an SNP Effect Matrix (SEM)
#' sem <- matrix(rnorm(12), ncol = 4)
#' colnames(sem) <- c("A", "C", "G", "T")
#' 
#' # create a list of SNPEffectMatrix objects
#' s <- SNPEffectMatrix(sem, baseline = -1, semId = "sem_id")
#' 
#' # plot the motif
#' plotSEM(s, semId(s))
#' 
#' # color by sequence
#' plotSEM(s, semId(s), motifSeq = "ACT", highlight = 2)
plotSEM <- function(s, motif=NULL, 
                    motifSeq=NULL, highlight=NULL, hcol="dodgerblue") {
  alt <- ref <- NA
  
  if ("SNPEffectMatrix" %in% is(s)) {
    sem_baseline <- baseline(s)
    sem_mtx <- sem(s)
    motif <- semId(s)
  } else if ("SNPEffectMatrixCollection" %in% is(s)) {
    # motif required if providing a SNPEffectMatrixCollection
    if (is.null(motif)) {
      stop("must specify the motif to plot if supplying a ",
      "SNPEffectMatrixCollection")
    }
    
    # check that motif given is in collection
    if (!(motif %in% names(sems(s)))) {
      stop("provided motif not found in SNPEffectMatrixCollection: ",
           "!(motif %in% names(sems(s))")
    }
    sem_baseline <- baseline(sems(s, semId = motif))
    sem_mtx <- sem(sems(s, semId = motif))
  } else {
    stop("s must be a SNPEffectMatrix or a SNPEffectMatrixCollection.")
  }
  
  # pivoting data longer
  bp <- rep(c("A", "C", "G", "T"), each=nrow(sem_mtx))
  motif_pos <- rep(1:nrow(sem_mtx), 4)
  sem_score <- sem_mtx |> 
    unlist()
  
  sem_mtx_long <- data.table(sem_score, motif_pos, bp)
  
  # if sequence is provided for plotting, 
  if (!is.null(motifSeq)) {
    if(nchar(motifSeq) != nrow(sem_mtx)) {
      stop(paste0("Number of characters in motifSeq (", nchar(motifSeq), 
                  ") does not equal the number of rows in the SEM (",
                  nrow(sem_mtx), ")"))
    }
    mseq <- unlist(strsplit(motifSeq, ""))
    
    mseq_dt <- data.table(mseq, 1:nrow(sem_mtx), mseq, "white") |> 
      stats::setNames(c("bp", "motif_pos", "mseq", "mcol"))
    
    sem_mtx_long <- merge(sem_mtx_long, mseq_dt, by=c("bp", "motif_pos"), 
                          all.x=T)
   
    sem_mtx_long$mcol[is.na(sem_mtx_long$mcol)] <- "#d7dbdd"
    sem_mtx_long$mseq[is.na(sem_mtx_long$mseq)] <- ""
    
  } 
  
  # if aren't plotting specific sequence
  if(is.null(motifSeq)) {
    sem_mtx_long$mcol <- "black"
  }
  
  text_size <- 7
  
  motif_plot <- ggplot2::ggplot(sem_mtx_long, 
                                aes(x=motif_pos, y=sem_score)) +
    theme_classic()
  
  if (!is.null(highlight)) {
    motif_plot <- motif_plot + geom_vline(xintercept = highlight, size = 15, 
                                          col=hcol, alpha = 0.1)
  }
  
  motif_plot <- motif_plot + 
    geom_hline(yintercept=0, color = "#eaeaea", size = 1) +
    geom_hline(yintercept=sem_baseline, linetype="dashed", 
               color = "#eaeaea", size = 1) + 
    geom_text(aes(label=sem_mtx_long$bp), size = text_size, col=sem_mtx_long$mcol) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = nrow(sem_mtx))) +
    theme(legend.position="none") +
    ggplot2::xlab("Position in Motif") +
    ggplot2::ylab("SNP Effect Score") +
    theme(text = element_text(size = 12)) +
    ggtitle(motif)
  
  
  if (!is.null(motifSeq)) {
    motif_plot <- motif_plot +
      geom_text(aes(label=mseq), size = text_size, 
                col=hcol, fontface="bold")
  }

  return(motif_plot)
}

