.semToPpm <- \(s) {
  norm_score <- apply(sem(s), 1, 
                      function(x) (2^x - 2^baseline(s)) / abs(2^baseline(s)))
  # replace negative scores with zero
  norm_score[norm_score < 0] <- 0
  ppm <- apply(norm_score, 2, function(x) x / sum(x))
  return(ppm)
}


.scrambleFlanks <- function(seqs) {
  prefix <- substr(seqs, 1, nchar(seqs)/2 - 0.5)
  suffix <- substr(seqs, nchar(prefix)+2, nchar(seqs))
  allele <- substr(seqs, nchar(prefix)+1, nchar(prefix)+1)
  
  flank_concat <- lapply(1:length(seqs), 
                         function(i) stringi::stri_c(prefix[i], suffix[i]))
  scrambled <- lapply(flank_concat, 
                      function(x) stringi::stri_rand_shuffle(x)) |> unlist()
  new_prefix <- substr(scrambled, 1, nchar(prefix))
  new_suffix <- substr(scrambled, nchar(prefix)+1, nchar(seqs)-1)
  
  scrambled_seqs <- stringi::stri_c(new_prefix, allele, new_suffix)
  return(scrambled_seqs)
}


#' Calculates binding enrichment for motif(s)
#' 
#' Perform a binomial test to determine if SNP Effect Matrices are bound more
#' often than expected.
#'
#' @param x `GRanges` object
#' @param x A `SNPEffectMatrix` or `SNPEffectMatrixCollection`
#' @param semList A `SNPEffectMatrix` or `SNPEffectMatrixCollection` object
#' @param bs_genome_obj A `BSgenome` object for the genome build to use.
#' @param background A list of DNA sequences to use for background. The length of each 
#' sequence must match the length of sequences in `x`
#'
#' @return a `list` of matrices
#'
#' @export
enrichSEMs <- \(x, semList, bs_genome_obj, background = NULL) {
  sem_names <- sems(semList) |> names()
  result <- data.frame(matrix(ncol = 3, nrow = length(sem_names)))
  colnames(result) <- c("pvalue", "n_successes", "bg_successes")
  rownames(result) <- sem_names
  
  sb <- scoreBinding(x, 
                     semList, 
                     bs_genome_obj, 
                     allele = "allele")
  if (is.null(background)) {
    offset <- lapply(sems(semList), function(x) {nrow(sem(x))}) |>
      unlist() |>
      max()
    x <- getFlankingSeqs(x = x, up = offset, down = offset, 
                          bs_genome_obj = bs_genome_obj, 
                          allele = "allele")
    
    sf <- .scrambleFlanks(x$seq)
    bg <- scoreBinding(sf, 
                       semList, 
                       BSgenome.Hsapiens.UCSC.hg38::Hsapiens)
  } else {
    bg <- scoreBinding(background, 
                       semList, 
                       BSgenome.Hsapiens.UCSC.hg38::Hsapiens, 
                       allele = "allele")
  }
  
  for (sem_name in sem_names) {
    successes <- sum(sb$scoreNorm[sb$sem_mtx_id == sem_name] > 0)
    sample_size <- length(sb$scoreNorm[sb$sem_mtx_id == sem_name])
    prob_success <- (sum(bg$scoreNorm[bg$sem_mtx_id == sem_name] > 0) + 1) / 
      length(bg$scoreNorm[bg$sem_mtx_id == sem_name])
    
    binom_test_res <- stats::binom.test(successes, sample_size, 
                                 p = prob_success, alternative = "greater")
    result[sem_name, "pvalue"] <- binom_test_res$p.value
    result[sem_name, "n_successes"] <- successes
    result[sem_name, "bg_successes"] <- sum(bg$scoreNorm[bg$sem_mtx_id == sem_name] > 0)
  }
  result$padj <- stats::p.adjust(result$pvalue, method = "BH")
  return(result)
}
