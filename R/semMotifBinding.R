#' Define a function to calculate risk/non-risk binding propensity
#' for all SEM motifs provided
#'
#' @param vr `VRanges` object with a single variant
#' @param sem_mtx a numeric matrix of motif position by nucleic acid
#' @param offset maximum number of offset basepairs offset from variant
#' @param baseline a numeric value of the baseline to use for normalization
scoreSem <- \(vr, sem_mtx, offset, baseline) {
  
  ## Score risk and non-risk sequences against each sem
  starts <- (offset-nrow(sem_mtx)+2):(offset+1)
  
  ref_seq <- vr$ref_seq
  alt_seq <- vr$alt_seq
  
  sem_mtx_t <- t(sem_mtx)
  
  nonRiskScores <-
    Biostrings::PWMscoreStartingAt(pwm = sem_mtx_t,
                                   subject = ref_seq,
                                   starting.at = starts)
  
  riskScores <-
    Biostrings::PWMscoreStartingAt(pwm = sem_mtx_t,
                                   subject = alt_seq,
                                   starting.at = starts)
  
  ## Find which frame has the best binding score (and record this)
  nonRiskFrame <- which.max(nonRiskScores)
  riskFrame <- which.max(riskScores)

  nonRiskScore <- nonRiskScores[nonRiskFrame]
  riskScore <- riskScores[riskFrame]

  ref_seq_frame <- substr(ref_seq,
                          start = starts[nonRiskFrame],
                          stop = starts[nonRiskFrame]+nrow(sem_mtx))

  alt_seq_frame <- substr(alt_seq,
                          start = starts[riskFrame],
                          stop = starts[riskFrame]+nrow(sem_mtx))


  nonRiskNorm <- (2^nonRiskScore - baseline) / abs(baseline)
  riskNorm <- (2^riskScore - baseline) / abs(baseline)
  
  scores <- data.table::data.table(nonRiskSeq=ref_seq_frame,
                                   riskSeq=alt_seq_frame,
                                   nonRiskScore=nonRiskScore,
                                   riskScore=riskScore,
                                   nonRiskNorm=nonRiskNorm,
                                   riskNorm=riskNorm)
  return(scores)
}

scoreVariant <- function(vr, semList) {
  varScore <- lapply(semList,
                     function(x) scoreSem(vr,
                                          sem_matrix(x),
                                          offset,
                                          baseline(x)))
  
  varScoreDt <- data.table::rbindlist(varScore)
  return(varScoreDt)
}


#' Define a function to calculate risk/non-risk binding propensity
#' for all SEM motifs provided
#'
#' @param vr `VRanges` object
#' @param semList a list of `SNPEffectMatrix` objects
#'
#' @return a SemplScores object
#'
#' @export
semMotifBinding <- \(vr, semList) {
  riskNorm <- riskSeq <- nonRiskNorm <- nonRiskSeq <- NULL
  
  ## Get maximum kmer length of all TFs ##
  offset <- max(unlist(lapply(semList, function(x) {nrow(sem_matrix(x))})))

  ## Collect up/downstream sequences
  vr <- get_flanking_seqs(vr, offset, offset)
   
  ## Create a new SemplScores object to store results
  semScores <- SemplScores(vr, sems = semList)
  
  new_scores <- lapply(seq_along(vr), 
                       function(i) scoreVariant(vr[i], semList)) |>
    rbindlist()
  
  scores(semScores)[, c("nonRiskSeq", "riskSeq", 
                        "nonRiskScore", "riskScore", 
                        "nonRiskNorm", "riskNorm")] <- new_scores

  return(semScores)
}

