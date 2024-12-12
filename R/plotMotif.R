#' Plot SEM scores for each nucleic acid in each position of the motif
#'
#' @param sems list of SNPEffectMatrix objects
#' @param motif sem id to plot. If providing a list of SNPEffectMatrix objects,
#' must match a name in the list. If providing a single SNPEffectMatrix object,
#' can provide any character string.
#' @param refSeq (optional) character sequence of reference
#' @param altSeq (optional) character sequence of alternate
#' @param refIndex (optional) index of variant location on reference
#' @param altIndex (optional) index of variant location on alternate
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
#' plotMotif(s, semId(s), refSeq = "ACT", altSeq = "ATT")
#' 
#' # color by sequence and highlight variant position
#' plotMotif(s, semId(s), refSeq = "ACT", altSeq = "ATT", refIndex = 2)
plotMotif <- function(sems, motif, 
                      refSeq=NULL, altSeq=NULL, 
                      refIndex=NULL, altIndex=NULL) {
  alt <- ref <- NA
  
  if (!is.null(refIndex) & !is.null(altIndex)) {
    if(refIndex != altIndex) {
      warning("ref and alt variant are located at different motif positions. \
              Suggest visualizing motifs seperately by sequence.")
    }
  }
  
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
  pwm_score <- sem_mtx |> 
    unlist()
  
  sem_mtx_long <- data.table(pwm_score, motif_pos, bp)
  
  if (!is.null(refSeq)) {
    text_color <- "grey"
    
    refseq <- unlist(strsplit(refSeq, ""))
    
    refseq_dt <- data.table(refseq, 1:nrow(sem_mtx), refseq) |> 
      setNames(c("bp", "motif_pos", "ref"))
    
    sem_mtx_long <- merge(sem_mtx_long, refseq_dt, by=c("bp", "motif_pos"), all.x=T)
    sem_mtx_long[is.na(sem_mtx_long)] <- ""
  } 
  
  if(!is.null(altSeq)) {    
    text_color <- "grey"

    altseq <- unlist(strsplit(altSeq, ""))
    
    altseq_dt <- data.table(altseq, 1:nrow(sem_mtx), altseq) |> 
      setNames(c("bp", "motif_pos", "alt"))
    
    sem_mtx_long <- merge(sem_mtx_long, altseq_dt, by=c("bp", "motif_pos"), all.x=T)
    sem_mtx_long[is.na(sem_mtx_long)] <- ""
    
  }
  
  if(is.null(refSeq) && is.null(altSeq)) {
    text_color <- "black"
  }
  
  text_size <- 8
  
  motif_plot <- ggplot2::ggplot(sem_mtx_long, 
                                aes(x=motif_pos, y=pwm_score))
  
  if (!is.null(refIndex)) {
    motif_plot <- motif_plot + geom_vline(xintercept = refIndex, size = 15, 
                                          col="blue", alpha = 0.1)
  }
  
  if (!is.null(altIndex)) {
    motif_plot <- motif_plot + geom_vline(xintercept = altIndex, size = 15, 
                                          col="red", alpha = 0.1)
  }
  
  motif_plot <- motif_plot + 
    geom_hline(yintercept=0, color = "#eaeaea", size = 1) +
    geom_hline(yintercept=sem_baseline, linetype="dashed", 
               color = "#eaeaea", size = 1) + 
    geom_text(aes(label=bp), size = text_size, col=text_color) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = nrow(sem_mtx))) +
    theme_classic() +
    theme(legend.position="none") +
    xlab("Position in Motif") +
    ylab("SNP Effect Score") +
    theme(text = element_text(size = 12)) +
    ggtitle(motif)
  
  if (!is.null(altSeq)) {
    motif_plot <- motif_plot +
      geom_text(aes(label=alt), size = text_size, 
                col="red", alpha=0.6, fontface="bold")
  }
  
  if (!is.null(refSeq)) {
    motif_plot <- motif_plot +
      geom_text(aes(label=ref), size = text_size, 
                col="blue", alpha=0.5, fontface="bold")
  }

  return(motif_plot)
}

