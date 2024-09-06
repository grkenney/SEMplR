#' Define a function to calculate risk/non-risk binding propensity
#' for all SEM motifs provided
#'
#' @param chr character of chromosome
#' @param pos integer of SNP position
#' @param ref character of reference nucleotide (can be upper or lower)
#' @param alt character of alternate nulceotide (can be upper or lower)
#' @param semFiles path to sem files
#' @param semBaselines path to sem baselines file
#'
semMotifBinding <- \(vr,
                     semFiles,
                     semBaselines) {

  ## Load SEM files ##
  ## Load libraries
  library(BSgenome)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(stringr)
  source(here::here("R/loadSEMs.R"))

  sems <- loadSEMs(semFiles)
  bl <- read.table(semBaselines)
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
        PWMscoreStartingAt(pwm = t(sems[[j]]),
                           subject = nonRiskSeq,
                           starting.at = (offset-nrow(sems[[j]])+2):(offset+1))
      riskScores <-
        PWMscoreStartingAt(pwm = t(sems[[j]]),
                           subject = riskSeq,
                           starting.at = (offset-nrow(sems[[j]])+2):(offset+1))

      ## Find which frame has the best binding score (and record this)
      if (max(nonRiskScores) > max(riskScores)){
        frame <- which.max(nonRiskScores)
      } else {
        frame <- which.max(riskScores)
      }

      semScores@scores <- rbind(semScores@scores,
                                data.table(seqnames=as.character(seqnames(semScores@variants[i])),
                                           ranges=start(semScores@variants[i]),
                                           sem=names(sems[j]),
                                           nonRiskSeq=str_sub(as.character(nonRiskSeq),
                                                              start = frame,
                                                              end = frame+nrow(sems[[j]])),
                                           riskSeq=str_sub(as.character(riskSeq),
                                                           start = frame,
                                                           end = frame+nrow(sems[[j]])),
                                           nonRiskScore=nonRiskScores[frame],
                                           riskScore=riskScores[frame],
                                           nonRiskNorm=(2^nonRiskScores[frame] - threshold[j]) / abs(threshold[j]),
                                           riskNorm=(2^riskScores[frame] - threshold[j]) / abs(threshold[j])))
    }
  }

  ## Return
  return(dt)
}
