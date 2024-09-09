#' Define a function to calculate risk/non-risk binding propensity
#' for all SEM motifs provided
#'
#' @param vr `VRanges` object
#' @param semFiles path to sem files
#' @param semBaselines path to sem baselines file
#'
#' @export
semMotifBinding <- \(vr,
                     semFiles,
                     semBaselines) {

  riskNorm <- riskSeq <- nonRiskNorm <- nonRiskSeq <- NULL

  ## Load SEM files ##
  sems <- loadSEMs(semFiles)
  bl <- utils::read.table(semBaselines)
  threshold <- 2^bl[match(gsub(".sem", "", basename(semFiles)), bl[,1]),2]

  ## Get maximum kmer length of all TFs ##
  offset <- max(unlist(lapply(sems, nrow)))

  ## Collect up/downstream sequences
  vr <- get_flanking_seqs(vr, offset, offset)

  ## Create a new SEMplR object to store results
  semScores <- SemplR(vr)

  for (i in 1:length(vr)){
    for(j in 1:length(sems)) {

      ## Score risk and non-risk sequences against each sem
      nonRiskScores <-
        Biostrings::PWMscoreStartingAt(pwm = t(sems[[j]]),
                           subject = vr[i]$ref_seq,
                           starting.at = (offset-nrow(sems[[j]])+2):(offset+1))
      riskScores <-
        Biostrings::PWMscoreStartingAt(pwm = t(sems[[j]]),
                           subject = vr[i]$alt_seq,
                           starting.at = (offset-nrow(sems[[j]])+2):(offset+1))

      ## Find which frame has the best binding score (and record this)
      if (max(nonRiskScores) > max(riskScores)){
        frame <- which.max(nonRiskScores)
      } else {
        frame <- which.max(riskScores)
      }

      semScores@scores <- rbind(semScores@scores,
                                data.table::data.table(seqnames=as.character(GenomeInfoDb::seqnames(semScores@variants[i])),
                                           ranges=BiocGenerics::start(semScores@variants[i]),
                                           sem=names(sems[j]),
                                           nonRiskSeq=stringr::str_sub(vr[i]$ref_seq,
                                                              start = frame,
                                                              end = frame+nrow(sems[[j]])),
                                           riskSeq=stringr::str_sub(vr[i]$alt_seq,
                                                           start = frame,
                                                           end = frame+nrow(sems[[j]])),
                                           nonRiskScore=nonRiskScores[frame],
                                           riskScore=riskScores[frame],
                                           nonRiskNorm=(2^nonRiskScores[frame] - threshold[j]) / abs(threshold[j]),
                                           riskNorm=(2^riskScores[frame] - threshold[j]) / abs(threshold[j])))
    }
  }

  ## Return
  return(semScores)
}
