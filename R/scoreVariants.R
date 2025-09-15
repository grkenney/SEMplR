.scoreAllele <- \(x, sem, prefix, alleleCol, nFlank, id) {
  score_col_suffixes <- c("Score", "Norm", "VarIndex", "Seq")
  score_cols <- paste0(prefix, score_col_suffixes)
  ds <- S4Vectors::mcols(x[, alleleCol]) |> unlist() |> unname()
  
  s <- lapply(sems(sem), 
              function(y) scoreSequence(sem = as.matrix(getSEM(y)), 
                                        dna_sequences = ds, 
                                        nFlank = nFlank, 
                                        bl = getBaseline(y),
                                        seqIds = id)) |>
    data.table::rbindlist(idcol = "SEM") |>
    stats::setNames(c("semId", score_cols, "varId"))
  return(s)
}

#' Calculate risk/non-risk binding propensity for all SEM motifs and
#' variants provided
#'
#' @param x `VRanges` object
#' @param sem a list of `SNPEffectMatrix` objects
#' @param genome A `BSgenome` object for the genome build to use. ie.
#' `BSgenome.Hsapiens.UCSC.hg19::Hsapiens`
#' @param refCol If providing a GRanges, the meta data column name with the 
#' reference (ref) allele. Ignored if providing a VRanges object.
#' @param altCol If providing a GRanges, the meta data column name with the
#' alternative (alt) allele. Ignored if providing a VRanges object.
#' @param varId A column name in the meta data of x to use as a unique id.
#'
#' @return a SEMplScores object
#'
#' @export
#' 
#' @examples
#' library(VariantAnnotation)
#'
#' # load default SEMs
#' data(sc)
#' 
#' # create a VRanges object
#' x <- VRanges(seqnames = "chr12",
#'               ranges = 94136009, 
#'               ref = "G", alt = "C")
#' 
#' # calculate binding propensity
#' scoreVariants(x, sc, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
#' 
scoreVariants <- \(x, sem, genome, 
                   refCol = NULL, altCol = NULL, 
                   varId = NULL) {
  riskNorm <- riskSeq <- nonRiskNorm <- nonRiskSeq <- NULL
  
  # Convert sem to a collection if it isn't one already
  sem <- .convertToSNPEffectMatrixCollection(sem)
  
  # Get maximum kmer length of all TFs ##
  nFlank <- lapply(sems(sem), function(x) {nrow(getSEM(x))}) |>
    unlist() |>
    max()

  # Collect up/downstream sequences
  x <- getRangeSeqs(x, genome = genome, 
                    up = nFlank, down = nFlank, 
                    refCol = refCol, altCol = altCol)
  
  # If a unique id
  if (is.null(varId)) {
    x$id <- lapply(seq_along(x), 
                   function(i) .makeVariantId(x[i], 
                                              refCol = refCol, 
                                              altCol = altCol)) |>
      unlist()
    id <- x$id
  } else {
    id <- S4Vectors::mcols(x)[, varId]
  }

  # Score each allele
  ref_scores <- .scoreAllele(x = x, sem = sem, 
                             prefix = "ref", alleleCol = "ref_seq", 
                             nFlank = nFlank, id = id)
  alt_scores <- .scoreAllele(x = x, sem = sem, 
                             prefix = "alt", alleleCol = "alt_seq", 
                             nFlank = nFlank, id = id)
  
  scores_merge <- merge(ref_scores, alt_scores, by = c("varId", "semId"))
  
  # reorder columns
  scores_merge <- scores_merge[, c("varId", "semId",
                                   "refSeq", "altSeq",
                                   "refScore", "altScore",
                                   "refNorm", "altNorm",
                                   "refVarIndex", "altVarIndex")]
  data.table::setkey(scores_merge, NULL) # clear the merge keys

  ## Store results in a SEMplScores object
  ss <- SEMplScores(ranges = x, 
                    semData = semData(sem), 
                    scores = scores_merge)

  return(ss)
}
