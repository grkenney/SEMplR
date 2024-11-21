#' Plot motif SEM scores
#'
#' @param sems list of SNPEffectMatrix objects
#' @param motif sem id to plot
#'
#' @import ggplot2
#'
#' @return a `ggplot` with sem scores for each amino acid per position
#'
#' @export
plotMotif <- function(sems, motif) {
  sem_baseline <- baseline(sems[[motif]])
  
  sem_mtx <- sem(sems[[motif]])
  sem_mtx$pos <- rownames(sem_mtx)
  
  aa_text_size <- 6
  
  motif_plot <- ggplot2::ggplot(sem_mtx, aes(x=pos, y=A)) +
    geom_hline(yintercept=0, color = "grey", size = 1) +
    geom_hline(yintercept=sem_baseline, linetype="dashed", 
               color = "grey", size = 1) + 
    geom_text(aes(label="A"), size = aa_text_size, col="#D81B60") +
    geom_text(aes(y = T, label="T"), size = aa_text_size, col="#1E88E5") +
    geom_text(aes(y = C, label="C"), size = aa_text_size, col="#004D40") +
    geom_text(aes(y = G, label="G"), size = aa_text_size, col="#A189D2") +
    theme_classic() +
    theme(legend.position="none") +
    xlab("Position in Motif") +
    ylab("SNP Effect Score") +
    theme(text = element_text(size = 12))
  return(motif_plot)
}
