
# get the start position of each sequence in concatenated string
getSeqStarts <- function(seq_list) {
  if (length(seq_list) == 1) {
    seq_starts <- 1
  } else {
    seq_lengths <- nchar(seq_list)[1:(length(seq_list)-1)]
    seq_starts <- append(1, seq_lengths) |> cumsum()
  }
  return(seq_starts)
}


# calculate the starting position of each frame to score in a single sequence
calcFrameStarts <- function(nbp, nflank, motif_size) {
  first_start <- nflank + 1 - (motif_size-1)
  last_start <- nflank + (nbp - 2*nflank)
  first_start:last_start
  return(first_start:last_start)
}


# concatenate sequences
concatSeqs <- function(seq_list) {
  # concatenate all variant sequences together
  concat_seq <- paste(seq_list, collapse="")
  return(concat_seq)
}


# adjust frame starting positions to reflect position in concatenated string
adjustFrameStarts <- function(frame_starts_list, seq_starts) {
  adjusted_frame_starts <- lapply(1:length(seq_starts), 
                                  function(i) 
                                    frame_starts_list[[i]] + seq_starts[i] - 1)
  return(adjusted_frame_starts)
}


fillNestedList <- function(fs, pwm_scores) {
  re <- lapply(fs, length) |> 
    unlist() |> 
    cumsum()
  
  rs <- lapply(1:length(fs), 
                         function(i) re[i]-length(fs[[i]])+1)
  nested_scores <- lapply(1:length(fs), 
                          function(i) pwm_scores[rs[[i]]:re[[i]]])
  return(nested_scores)
}


# score multiple frames within a sequence against an sem matrix
scoreMatrix <- function(sem_mtx, concat_seq, frame_starts) {
  # score all sequences and all frames against the pwm
  pwm_scores <- Biostrings::PWMscoreStartingAt(pwm = t(sem_mtx), 
                                   subject = concat_seq, 
                                   starting.at = unlist(frame_starts))
  nested_pwm_scores <- fillNestedList(frame_starts, pwm_scores)
  return(nested_pwm_scores)
}


# normalize scoring against baselines
normMatrix <- function(max_scores, bl) {
  score_norm <- (2^max_scores - 2^bl) / abs(2^bl)
  return(score_norm)
}


#' Define a function to calculate risk/non-risk binding propensity
#' for all SEM motifs provided
#'
#' @param varId a list of unique character identifiers for each variant
#' @param varSeq a list of variant sequences
#' @param semObj SNPEffectMatrix object
#' @param nflank maximum number of offset basepairs offset from variant
scoreVariants <- function(varId, varSeq, semObj, nflank) {
  
  sm <- sem(semObj)
  bl <- baseline(semObj)
  
  # count number of variants and number of frames (number of sem positions)
  nVars <- length(varSeq)
  frame_len <- nrow(sm)
  
  if (min(nchar(varSeq)) < frame_len) {
    stop("Variant sequence length ", min(nchar(varSeq)), 
         " is less than SEM length of ", frame_len, 
         ". \n Make sure all variant sequences are greater than or equal",
         "to SEM length")
  }
  
  frame_starts <- lapply(varSeq, 
                         function(x) 
                           calcFrameStarts(nbp = nchar(x), 
                                           nflank = nflank, 
                                           motif_size = frame_len))
  
  cseq <- concatSeqs(varSeq)
  seq_starts <- getSeqStarts(seq_list = varSeq)
  adjusted_frames <- adjustFrameStarts(frame_starts, seq_starts)
  
  scoresList <- scoreMatrix(sm, cseq, adjusted_frames)
  
  maxScoreIndex <- lapply(scoresList, which.max)
  maxScores <- lapply(1:length(scoresList), 
                      function(i) scoresList[[i]][[maxScoreIndex[[i]]]]) |>
    unlist()
  
  # get the top scoring frame
  frameSeq <- lapply(1:length(adjusted_frames), 
                     function(i) substr(cseq, 
                                        adjusted_frames[[i]][[maxScoreIndex[[i]]]], 
                                        adjusted_frames[[i]][[maxScoreIndex[[i]]]] + frame_len - 1)) |>
    unlist()
  
  maxFrameStart <- lapply(1:length(frame_starts),
                          function(i) frame_starts[[i]][maxScoreIndex[[i]]]) |>
    unlist()
  
  normScores <- normMatrix(maxScores, bl)
  
  scores_dt <- data.table::data.table(var_id=varId,
                                      sem_mtx_id=semId(semObj),
                                      seq=frameSeq,
                                      vari=maxFrameStart,
                                      score=maxScores,
                                      scoreNorm=normScores)
  
  return(scores_dt)
}


#' Calculate risk/non-risk binding propensity for all SEM motifs and
#' variants provided
#'
#' @param vr `VRanges` object
#' @param semList a list of `SNPEffectMatrix` objects
#' @param bs_genome_obj A `BSgenome` object for the genome build to use.
#'
#' @return a SemplScores object
#'
#' @export
#' 
#' @examples
#' library(VariantAnnotation)
#'
#' # create an SNP Effect Matrix (SEM)
#' sem <- matrix(rnorm(12), ncol = 4)
#' colnames(sem) <- c("A", "C", "G", "T")
#' 
#' # create a VRanges object
#' vr <- VRanges(seqnames = "chr12",
#'               ranges = 94136009, 
#'               ref = "G", alt = "C")
#' 
#' # create a list of SNPEffectMatrix objects
#' semList <- SNPEffectMatrix(sem, baseline = -1, semId = "sem_id")
#' 
#' # calculate binding propensity
#' semMotifBinding(vr, semList)
#' 
semMotifBinding <- \(vr, semList, 
                     bs_genome_obj=BSgenome.Hsapiens.UCSC.hg19::Hsapiens) {
  riskNorm <- riskSeq <- nonRiskNorm <- nonRiskSeq <- NULL
  
  # if semList is not a list, make it a named list
  if (!is.list(semList)) {
    semList <- list(semList) |> stats::setNames(semId(semList))
  }
  
  ## Get maximum kmer length of all TFs ##
  offset <- max(unlist(lapply(semList, function(x) {nrow(sem(x))})))

  ## Collect up/downstream sequences
  vr <- getFlankingSeqs(vr, offset, offset, bs_genome_obj)
   
  ## Create a new SemplScores object to store results
  semScores <- SemplScores(vr, semList)

  ref_scores <- lapply(semList,
                       function(x) scoreVariants(variants(semScores)$id, 
                                                 vr$ref_seq, x, offset) ) |>
    data.table::rbindlist()
  ref_scores_merge <- merge(ref_scores, scores(semScores), 
                            by.x = c("var_id", "sem_mtx_id"), 
                            by.y = c("varId", "semId"))
  scores(semScores)[, c("nonRiskSeq", "nonRiskVarIndex", "nonRiskScore", "nonRiskNorm")] <- ref_scores_merge[, c("seq", "vari", "score", "scoreNorm")]

  alt_scores <- lapply(semList,
                       function(x) scoreVariants(variants(semScores)$id,
                                                 vr$alt_seq, x, offset) ) |>
    data.table::rbindlist()
  
  alt_scores_merge <- merge(alt_scores, scores(semScores), 
                            by.x = c("var_id", "sem_mtx_id"), 
                            by.y = c("varId", "semId"))
  scores(semScores)[, c("riskSeq", "riskVarIndex", "riskScore", "riskNorm")] <- alt_scores_merge[, c("seq", "vari", "score", "scoreNorm")]

  return(semScores)
}

