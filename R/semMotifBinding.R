#' Define a function to calculate risk/non-risk binding propensity
#' for all SEM motifs provided
#'
#' @param vr `VRanges` object with a single variant
#' @param sem_matrix a numeric matrix of motif position by nucleic acid
#' @param tf_name a string of the transcription factor name
#' @param threshold a numeric value of the threshold to use for normalization
scoreVariant <- \(vr, sem_matrix, tf_name, offset, threshold) {

  ## Score risk and non-risk sequences against each sem
  starting_indices <- (offset-nrow(sem_matrix)+2):(offset+1)
  nonRiskScores <-
    Biostrings::PWMscoreStartingAt(pwm = t(sem_matrix),
                                   subject = vr$ref_seq,
                                   starting.at = starting_indices)
  riskScores <-
    Biostrings::PWMscoreStartingAt(pwm = t(sem_matrix),
                                   subject = vr$alt_seq,
                                   starting.at = starting_indices)

  ## Find which frame has the best binding score (and record this)
  nonRiskFrame <- which.max(nonRiskScores)
  riskFrame <- which.max(riskScores)

  ref_seq_frame <- stringr::str_sub(vr$ref_seq,
                                    start = starting_indices[nonRiskFrame],
                                    end = starting_indices[nonRiskFrame]+nrow(sem_matrix))

  alt_seq_frame <- stringr::str_sub(vr$alt_seq,
                                    start = starting_indices[riskFrame],
                                    end = starting_indices[riskFrame]+nrow(sem_matrix))

  nonRiskNorm <- (2^nonRiskScores[nonRiskFrame] - threshold) / abs(threshold)
  riskNorm <- (2^riskScores[riskFrame] - threshold) / abs(threshold)

  scores <- data.table::data.table(seqnames=as.character(GenomeInfoDb::seqnames(vr)),
                                   ranges=BiocGenerics::start(vr),
                                   sem=tf_name,
                                   nonRiskSeq=ref_seq_frame,
                                   riskSeq=alt_seq_frame,
                                   nonRiskScore=nonRiskScores[nonRiskFrame],
                                   riskScore=riskScores[riskFrame],
                                   nonRiskNorm=nonRiskNorm,
                                   riskNorm=riskNorm)
  return(scores)
}


#' Define a function to calculate risk/non-risk binding propensity
#' for all SEM motifs provided
#'
#' @param vr `VRanges` object
#' @param semList a list of `SNPEffectMatrix` objects. Defaults to loading .sems
#' and baseline data from github.com/Boyle-Lab/SEMpl/raw/master/SEMs/
#'
#' @export
semMotifBinding <- \(vr,
                     semList=NULL) {

  riskNorm <- riskSeq <- nonRiskNorm <- nonRiskSeq <- NULL

  ## Load SEM files ##
  if (is.null(semList)) {
    sems <- loadSEMs()
  } else {
    sems <- semList
  }

  threshold <- lapply(sems, function(x) {2^x@baseline})

  ## Get maximum kmer length of all TFs ##
  offset <- max(unlist(lapply(sems, function(x) {nrow(x@matrix)})))

  ## Collect up/downstream sequences
  vr <- get_flanking_seqs(vr, offset, offset)

  ## Create a new SEMplR object to store results
  semScores <- SemplR(vr)

  semScores@scores <- data.table::data.table()

  for (i in 1:length(vr)){
    varScore <- lapply(1:length(semList),
                       function(j) { scoreVariant(semScores@variants[i],
                                                  semList[[j]]@matrix,
                                                  semList[[j]]@tf_name,
                                                  offset,
                                                  threshold[[j]]) })
    varScoreDt <- data.table::rbindlist(varScore)

    semScores@scores <- rbind(semScores@scores, varScoreDt)
  }

  ## Return
  return(semScores)
}
