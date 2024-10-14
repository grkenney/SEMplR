#' Define a function to calculate risk/non-risk binding propensity
#' for all SEM motifs provided
#'
#' @param vr `VRanges` object with ref_seq and alt_seq metadata columns
#' @param sem_mtx a numeric matrix of motif position by nucleic acid
#' @param bl a numeric value of the baseline to use for normalization
#' @param offset maximum number of offset basepairs offset from variant
scoreVariants <- function(var_id, var_seq, sem_mtx_id, sem_mtx, bl, offset) {
  concat_seq <- paste(var_seq, collapse="")
  nvars <- length(var_seq)
  nframes <- nrow(sem_mtx)
  
  var_starts <- (offset*2+1) * (0:(nvars-1))
  
  frame_starts <- unlist( lapply( var_starts,
                                  function(n) n + (offset-nframes+2):(offset+1) ) )
  
  pwm_scores <- PWMscoreStartingAt(pwm = t(sem_mtx), 
                                   subject = concat_seq, 
                                   starting.at = frame_starts)
  scores_mtx <- matrix(pwm_scores, nrow = nframes)
  
  score_max <- colMaxs(scores_mtx)
  score_norm <- (2^score_max - 2^bl) / abs(2^bl)
  
  score_max_i <- apply(scores_mtx, 2, which.max)
  
  # get the sequence of the frame
  seq_frame_start <- unlist(lapply(1:length(score_max_i), 
                                   function(i) ( frame_starts[score_max_i[i] + (0:nvars * nframes)[i]] )))
  
  seq_frame <- unlist( lapply( seq_frame_start, 
                               function(x) substr(concat_seq, x, x + nframes - 1) ))
  
  scores_dt <- data.table::data.table(var_id=var_id,
                                      sem_mtx_id=sem_mtx_id,
                                      seq=seq_frame,
                                      score=score_max,
                                      scoreNorm=score_norm)
  
  return(scores_dt)
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
  
  ref_scores <- lapply(semList,
                       function(x) scoreVariants(vr$id, vr$ref_seq, sem_id(x), sem_matrix(x), baseline(x), offset) ) |>
    rbindlist()
  ref_scores_merge <- merge(ref_scores, scores(semScores), 
                            by.x = c("var_id", "sem_mtx_id"), 
                            by.y = c("varId", "sem"))
  scores(semScores)[, c("nonRiskSeq", "nonRiskScore", "nonRiskNorm")] <- ref_scores_merge[, c("seq", "score", "scoreNorm")]

  alt_scores <- lapply(semList,
                       function(x) scoreVariants(vr$id, vr$alt_seq, sem_id(x), sem_matrix(x), baseline(x), offset) ) |>
    rbindlist()
  
  alt_scores_merge <- merge(alt_scores, scores(semScores), 
                            by.x = c("var_id", "sem_mtx_id"), 
                            by.y = c("varId", "sem"))
  scores(semScores)[, c("riskSeq", "riskScore", "riskNorm")] <- alt_scores_merge[, c("seq", "score", "scoreNorm")]

  return(semScores)
}

