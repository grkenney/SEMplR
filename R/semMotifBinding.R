concatVars <- function(var_seq) {
  # concatenate all variant sequences together
  concat_seq <- paste(var_seq, collapse="")
  return(concat_seq)
}


calcVarStarts <- function(offset, nVars, nFrames) {
  # calculate the first position of each variant
  var_starts <- (offset*2+1) * (0:(nVars-1))
  
  # find all frame starts for all variants
  frame_starts <- unlist( lapply( var_starts,
                                  function(n) n + (offset-nFrames+2):(offset+1) ) )
  return(frame_starts)
}


scoreMatrix <- function(sem_mtx, concat_seq, frame_starts, nFrames) {
  # score all sequences and all frames against the pwm
  pwm_scores <- Biostrings::PWMscoreStartingAt(pwm = t(sem_mtx), 
                                   subject = concat_seq, 
                                   starting.at = frame_starts)
  
  scores_mtx <- matrix(pwm_scores, nrow = nFrames)
  return(scores_mtx)
}


normMatrix <- function(scores_mtx, bl) {
  score_max <- MatrixGenerics::colMaxs(scores_mtx)
  score_norm <- (2^score_max - 2^bl) / abs(2^bl)
  return(score_norm)
}


#' Define a function to calculate risk/non-risk binding propensity
#' for all SEM motifs provided
#'
#' @param varId a list of unique character identifiers for each variant
#' @param varSeq a list of variant sequences
#' @param semId a unique character identifier for the sem
#' @param semMtx a numeric matrix of motif position by nucleic acid
#' @param bl a numeric value of the baseline to use for normalization
#' @param offset maximum number of offset basepairs offset from variant
scoreVariants <- function(varId, varSeq, semId, semMtx, bl, offset) {
  # concatenate all variant sequences together
  concatSeq <- paste(varSeq, collapse="")
  
  # count number of variants and number of frames (number of sem positions)
  nVars <- length(varSeq)
  nFrames <- nrow(semMtx)
  
  # calculate position of frame starts in concatenated sequence
  frameStarts <- calcVarStarts(offset, nVars, nFrames)
  
  scoresMtx <- scoreMatrix(semMtx, concatSeq, frameStarts, nFrames)
  
  maxScoreIndex <- apply(scoresMtx, 2, which.max)
  
  # get the sequence of the frame
  seq_frame_start <- unlist(lapply(1:length(maxScoreIndex), 
                                   function(i) ( frameStarts[maxScoreIndex[i] + (0:nVars * nFrames)[i]] )))
  
  seq_frame <- unlist( lapply( seq_frame_start, 
                               function(x) substr(concatSeq, x, x + nFrames - 1) ))
  
  score_norm <- normMatrix(scoresMtx, bl)
  
  scores_dt <- data.table::data.table(var_id=varId,
                                      sem_mtx_id=semId,
                                      seq=seq_frame,
                                      score=scoresMtx[maxScoreIndex],
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
  offset <- max(unlist(lapply(semList, function(x) {nrow(sem(x))})))

  ## Collect up/downstream sequences
  vr <- get_flanking_seqs(vr, offset, offset)
   
  ## Create a new SemplScores object to store results
  semScores <- SemplScores(vr, sems = semList)
  
  ref_scores <- lapply(semList,
                       function(x) scoreVariants(vr$id, vr$ref_seq, semId(x), sem(x), baseline(x), offset) ) |>
    data.table::rbindlist()
  ref_scores_merge <- merge(ref_scores, scores(semScores), 
                            by.x = c("var_id", "sem_mtx_id"), 
                            by.y = c("varId", "semId"))
  scores(semScores)[, c("nonRiskSeq", "nonRiskScore", "nonRiskNorm")] <- ref_scores_merge[, c("seq", "score", "scoreNorm")]

  alt_scores <- lapply(semList,
                       function(x) scoreVariants(vr$id, vr$alt_seq, semId(x), sem(x), baseline(x), offset) ) |>
    data.table::rbindlist()
  
  alt_scores_merge <- merge(alt_scores, scores(semScores), 
                            by.x = c("var_id", "sem_mtx_id"), 
                            by.y = c("varId", "semId"))
  scores(semScores)[, c("riskSeq", "riskScore", "riskNorm")] <- alt_scores_merge[, c("seq", "score", "scoreNorm")]

  return(semScores)
}

