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
  }

  threshold <- lapply(sems, function(x) {2^x@baseline})

  ## Get maximum kmer length of all TFs ##
  offset <- max(unlist(lapply(sems, function(x) {nrow(x@matrix)})))

  ## Collect up/downstream sequences
  vr <- get_flanking_seqs(vr, offset, offset)

  ## Create a new SEMplR object to store results
  semScores <- SemplR(vr)

  for (i in 1:length(vr)){
    for(j in 1:length(sems)) {

      sem_matrix <- sems[[j]]@matrix

      ## Score risk and non-risk sequences against each sem
      nonRiskScores <-
        Biostrings::PWMscoreStartingAt(pwm = t(sem_matrix),
                           subject = vr[i]$ref_seq,
                           starting.at = (offset-nrow(sem_matrix)+2):(offset+1))
      riskScores <-
        Biostrings::PWMscoreStartingAt(pwm = t(sem_matrix),
                           subject = vr[i]$alt_seq,
                           starting.at = (offset-nrow(sem_matrix)+2):(offset+1))

      ## Find which frame has the best binding score (and record this)
      if (max(nonRiskScores) > max(riskScores)){
        frame <- which.max(nonRiskScores)
      } else {
        frame <- which.max(riskScores)
      }

      semScores@scores <- rbind(semScores@scores,
                                data.table::data.table(seqnames=as.character(GenomeInfoDb::seqnames(semScores@variants[i])),
                                           ranges=BiocGenerics::start(semScores@variants[i]),
                                           sem=sems[[j]]@tf_name,
                                           nonRiskSeq=stringr::str_sub(vr[i]$ref_seq,
                                                              start = frame,
                                                              end = frame+nrow(sem_matrix)),
                                           riskSeq=stringr::str_sub(vr[i]$alt_seq,
                                                           start = frame,
                                                           end = frame+nrow(sem_matrix)),
                                           nonRiskScore=nonRiskScores[frame],
                                           riskScore=riskScores[frame],
                                           nonRiskNorm=(2^nonRiskScores[frame] - threshold[[j]]) / abs(threshold[[j]]),
                                           riskNorm=(2^riskScores[frame] - threshold[[j]]) / abs(threshold[[j]])))
    }
  }

  ## Return
  return(semScores)
}
