# given a GRanges object, build a position identifier string of format
# seqname:position
.makePositionId <- \(x) {
  start_pos <- IRanges::start(IRanges::ranges(x))
  end_pos <- IRanges::end(IRanges::ranges(x))
  sn <- GenomeInfoDb::seqnames(x)
  pos_str <- ifelse(start_pos == end_pos, 
                    start_pos, paste0(start_pos, "-", end_pos))
  
  id <- paste0(sn, ":", pos_str)
  return(id)
}

# # score a single GRanges position against all SEMs
# .scoreSingleSEM <- \(x, semList, bs_genome_obj, offset) {
#   id <- .makePositionId(x)
#   s <- lapply(sems(semList), 
#               function(y) calculateScores(id, varSeq = x$seq, 
#                                           semObj = y, 
#                                           nflank = offset)) |>
#     data.table::rbindlist()
#   return(s)
# }


.test_if_sequence_list <- \(x) {
  is_sequence_list <- T
  # if is neither a vector or a list, it's not collection of DNA seqs
  if (!(is.vector(x)) & !(is.list(x))) {
    is_sequence_list <- F
  }
  # check if it's a character vector
  if (!is.character(x)) {
    is_sequence_list <- F
  }
  return(is_sequence_list)
}


#' Calculate binding propensity for all SEM motifs and
#' genomic positions provided
#'
#' @param x `GRanges` object
#' @param semList A `SNPEffectMatrix` or `SNPEffectMatrixCollection` object
#' @param bs_genome_obj A `BSgenome` object for the genome build to use. ie.
#' `BSgenome.Hsapiens.UCSC.hg19::Hsapiens`
#' @param allele Column name in meta data storing allele
#'
#' @return a `data.table` object
#'
#' @export
#' 
#' @examples
#' library(GenomicRanges)
#'
#' # create an SNP Effect Matrix (SEM)
#' s <- matrix(rnorm(12), ncol = 4)
#' colnames(s) <- c("A", "C", "G", "T")
#' 
#' # create a GRanges object
#' gr <- GenomicRanges::GRanges(seqnames = "chr12",
#'               ranges = 94136009, 
#'               ref = "G", alt = "C")
#' 
#' # create a list of SNPEffectMatrix objects
#' semList <- SNPEffectMatrix(s, baseline = -1, semId = "sem_id")
#' 
#' # calculate binding propensity
#' scoreBinding(gr, semList, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
#' 
scoreBinding <- \(x, semList, bs_genome_obj, allele = NULL) {
  is_sequence_list <- .test_if_sequence_list(x)
  if(is(x)[1] != "GRanges" & !is_sequence_list) {
    stop(is(x)[1], " is not an accepted class for scoreBinding. 
         x must be a GRanges object or vector of DNA sequences")
  }
  
  # if given a SNPEffectMatrix or a list of SNPEffectMatrix objects,
  # make it a SNPEffectMatrixCollection
  if ("SNPEffectMatrix" %in% is(semList)) {
    semList <- SNPEffectMatrixCollection(list(semList))
  } else if (is(semList)[1] == "list") {
    if (is(semList[[1]]) == "SNPEffectMatrix") {
      semList <- SNPEffectMatrixCollection(semList)
      }
  } else if (!("SNPEffectMatrixCollection" %in% is(semList))) {
    stop("semList must be a SNPEffectMatrixCollection object, 
         SNPEffectMatrix object, or a list of SNPEffectMatrix objects. 
         \nSee ?SNPEffectMatrixCollection or use the provided default `sc`")
  }
  
  ## Get maximum kmer length of all TFs ##
  offset <- lapply(sems(semList), function(x) {nrow(sem(x))}) |>
    unlist() |>
    max()
  
  # get flanking sequences
  if (!is_sequence_list) {
    x <- getFlankingSeqs(x, offset, offset, bs_genome_obj, allele)
    seqs <- x$seq
    id <- lapply(seq_along(x), function(i) .makePositionId(x = x[i])) |>
      unlist()
  } else {
    seqs <- x
    id <- 1:length(seqs)
  }
  
  s <- lapply(sems(semList), 
              function(y) calculateScores(id, varSeq = seqs, 
                                          semObj = y, 
                                          nflank = offset)) |>
    data.table::rbindlist()

  return(s)
}
