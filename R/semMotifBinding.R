#' Define a function to calculate risk/non-risk binding propensity
#' for all SEM motifs provided
#'
#' @param vr `VRanges` object with a single variant
#' @param sem_matrix a numeric matrix of motif position by nucleic acid
#' @param sem_id a string with the unique id of the SEM
#' @param threshold a numeric value of the threshold to use for normalization
scoreVariant <- \(vr, sem_matrix, sem_id, offset, threshold) {

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
                                   sem=sem_id,
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
#' @return a SemplScores object
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

  ## Create a new SemplScores object to store results
  semScores <- SemplScores(vr)
  semScores@sem_metadata <- lapply(semList, 
                                   function(x) data.table::data.table(sem_id=x@sem_id,
                                                          tf_name=x@tf_name,
                                                          ensembl_id=x@ensembl_id,
                                                          uniprot_id=x@uniprot_id,
                                                          cell_type=x@cell_type)) |>
    data.table::rbindlist()

  semScores@scores <- data.table::data.table()

  for (i in 1:length(vr)){
    varScore <- lapply(1:length(semList),
                       function(j) { scoreVariant(semScores@variants[i],
                                                  semList[[j]]@matrix,
                                                  semList[[j]]@sem_id,
                                                  offset,
                                                  threshold[[j]]) })
    varScoreDt <- data.table::rbindlist(varScore)

    semScores@scores <- rbind(semScores@scores, varScoreDt)
  }

  ## Return
  return(semScores)
}
