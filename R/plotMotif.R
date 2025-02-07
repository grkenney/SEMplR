#' Plot SEM scores for each nucleic acid in each position of the motif
#'
#' @param sems list of SNPEffectMatrix objects
#' @param motif sem id to plot. If providing a list of SNPEffectMatrix objects,
#' must match a name in the list. If providing a single SNPEffectMatrix object,
#' can provide any character string.
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
#' plotMotif(s, semId(s))
#' 
#' # color by sequence
#' plotMotif(s, semId(s), motifSeq = "ACT", highlight = 2)
plotMotif <- function(sems, motif, 
                      motifSeq=NULL, highlight=NULL, hcol="dodgerblue") {
  alt <- ref <- NA
  
  if (is.list(sems)) {
    sem_baseline <- baseline(sems[[motif]])
    sem_mtx <- sem(sems[[motif]])
  } else {
    sem_baseline <- baseline(sems)
    sem_mtx <- sem(sems)
  }
  
  # pivoting data longer
  bp <- rep(c("A", "C", "G", "T"), each=nrow(sem_mtx))
  motif_pos <- rep(1:nrow(sem_mtx), 4)
  sem_score <- sem_mtx |> 
    unlist()
  
  sem_mtx_long <- data.table(sem_score, motif_pos, bp)
  
  # if sequence is provided for plotting, 
  if (!is.null(motifSeq)) {
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
    xlab("Position in Motif") +
    ylab("SNP Effect Score") +
    theme(text = element_text(size = 12)) +
    ggtitle(motif)
  
  
  if (!is.null(motifSeq)) {
    motif_plot <- motif_plot +
      geom_text(aes(label=mseq), size = text_size, 
                col=hcol, fontface="bold")
  }

  return(motif_plot)
}

