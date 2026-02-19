.scoreAllele <- function(x, sem, prefix, alleleCol, nFlank, genome, id, rc) {
    score_col_suffixes <- c("Score", "Norm", "VarIndex", "Seq")
    score_cols <- paste0(prefix, score_col_suffixes)
    ds <- S4Vectors::mcols(x[, alleleCol]) |>
        unlist() |>
        unname()

    s <- scoreBinding(ds,
        sem = sem, genome = genome,
        nFlank = nFlank, seqId = id, rc = rc
    )
    colnames(s) <- c("varId", "SEM", "rc", score_cols)
    return(s)
}


.getVarId <- function(x, varId, refCol, altCol) {
    if (is.null(varId)) {
        x$id <- lapply(
            seq_along(x),
            function(i) {
                .makeVariantId(x[i],
                               refCol = refCol,
                               altCol = altCol
                )
            }
        ) |>
            unlist()
        id <- x$id
    } else {
        id <- S4Vectors::mcols(x)[, varId]
    }
    return(id)
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
#' @param rc plot the reverse complement SEMs
#'
#' @return a SEMScores object with slots for the ranges of the provided
#' variants with an added `sequence` column with the sequence scored, 
#' the SEM metadata, and the resulting scoring table.
#' These slots are accessible with the `getRanges()`, `semData()`, and 
#' `scores()` accessor functions respectively.
#' 
#' The scoring table will contain the following columns:
#' - varId: a unique identifier for the sequence scored
#' - SEM: the identifier for the SEM
#' - rc: the orientation of the sequence (fwd or rev)
#' - score: the unnormalized SEM score
#' - scoreNorm: the SEM score normalized to it's corresponding baseline
#' - index: the index of the motif with the highest SEM score within the 
#' sequence scored
#' - seq: the motif sequence with the highest SEM scores within the sequence
#' scored
#'
#' @export
#'
#' @examples
#' library(VariantAnnotation)
#'
#' # load default SEMs
#'
#' # create a VRanges object
#' x <- VRanges(
#'     seqnames = "chr12",
#'     ranges = 94136009,
#'     ref = "G", alt = "C"
#' )
#'
#' # calculate binding propensity
#' scoreVariants(x, SEMC, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
#'
scoreVariants <- function(x, sem, genome,
    refCol = NULL, altCol = NULL,
    varId = NULL, rc = TRUE) {
    riskNorm <- riskSeq <- nonRiskNorm <- nonRiskSeq <- NULL

    # Convert sem to a collection if it isn't one already
    sem <- .convertToSNPEffectMatrixCollection(sem)

    # Get maximum kmer length of all TFs ##
    nFlank <- lapply(getSEMs(sem), function(x) {
        nrow(getSEM(x))
    }) |> unlist() |> max()

    # Collect up/downstream sequences
    x <- getRangeSeqs(x,
        genome = genome,
        up = nFlank, down = nFlank,
        refCol = refCol, altCol = altCol
    )

    id <- .getVarId(x, varId, refCol, altCol)

    # Score each allele
    ref_scores <- .scoreAllele(
        x = x, sem = sem,
        prefix = "ref", alleleCol = "ref_seq",
        nFlank = nFlank, genome = genome, id = id, rc = rc
    )
    alt_scores <- .scoreAllele(
        x = x, sem = sem,
        prefix = "alt", alleleCol = "alt_seq",
        nFlank = nFlank, genome = genome, id = id, rc = rc
    )

    scores_merge <- merge(ref_scores, alt_scores, by = c("varId", "SEM", "rc"))

    # reorder columns
    scores_merge <- scores_merge[, c(
        "varId", "SEM", "rc", "refSeq", "altSeq",
        "refScore", "altScore", "refNorm", "altNorm",
        "refVarIndex", "altVarIndex"
    )]
    data.table::setkey(scores_merge, NULL) # clear the merge keys

    ## Store results in a SEMScores object
    ss <- SEMScores(ranges = x, semData = semData(sem), scores = scores_merge)

    return(ss)
}
