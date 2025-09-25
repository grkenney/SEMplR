#' Convert a SNP Effect Matrix to a position probability matrices (PPMs)
#'
#' Converts the SNP Effect Matrix in a SNPEffectMatrix or
#' SNPEffectMatrixCollection object to a position probability matrix (PPM)
#' where each entry is the probability of a base at each position in the motif.
#'
#' @param x A `SNPEffectMatrix` or `SNPEffectMatrixCollection` object
#'
#' @importFrom methods is
#'
#' @return a `list` of matrices
#'
#' @export
#'
#' @examples
#' # Load default SEMs
#' data(SEMC)
#'
#' # From a SNPEffectMatrixCollection
#' convertSEMsToPPMs(SEMC)
#'
#' # From a single SNPEffectMatrix
#' convertSEMsToPPMs(getSEMs(SEMC, "AP2B_HUMAN.SK-N-SH"))
#'
convertSEMsToPPMs <- \(x) {
    if (is(x, "SNPEffectMatrixCollection")) {
        ss <- getSEMs(x)
    } else if (is(x, "SNPEffectMatrix")[1]) {
        ss <- list(x)
    } else if (is(x, "list")[1] & is(x[[1]], "SNPEffectMatrix")) {
        ss <- x
    } else {
        rlang::abort(paste0(
            "x must be of class SNPEffectMatrixCollection or ",
            "SNPEffectMatrix or a list of SNPEffectMatrix objects"
        ))
    }

    ppms <- lapply(ss, function(s) .semToPpm(s))
    return(ppms)
}
