# ---- constructor ----

#' SEMScores object and constructor
#'
#' Constructs a SEMScores class object.
#'
#' @param ranges A `GRanges` or `VRanges` object to hold one or more variants
#' @param semData A named list of SNPEffectMatrix objects
#' @param scores (optional) A `data.table` object for motif information and
#' binding scores
#'
#' @importFrom methods new
#' @importFrom VariantAnnotation VRanges
#' @importFrom S4Vectors mcols
#'
#' @return a SEMScores object
#' @docType class
#' @export
#'
#' @examples
#' # load default SEMs
#' data(SEMC)
#'
#' # create a VRanges object
#' vr <- VariantAnnotation::VRanges(
#'     seqnames = c("chr12", "chr19"),
#'     ranges = c(94136009, 10640062),
#'     ref = c("G", "T"), alt = c("C", "A")
#' )
#'
#' SEMScores(ranges = vr, semData = semData(SEMC))
#'
SEMScores <- function(ranges = NULL, semData = NULL, scores = NULL) {
    # if no ranges provided, make an empty VRanges object
    if (all(is.null(ranges))) {
        r <- VariantAnnotation::VRanges()
    } else {
        r <- ranges
    }

    if (is.null(scores)) {
        scores_table <- data.table()
    } else {
        scores_table <- scores
    }

    if (is.null(semData)) {
        semData <- data.table()
    }

    new("SEMScores",
        ranges = r,
        semData = semData,
        scores = scores_table
    )
}


# ---- accessors ----

#' Access ranges slot in a SEMScores object
#'
#' @param x a SEMScores object
#' @rdname getRanges
#' @export
#'
#' @return A GRanges or VRanges object
#'
#' @examples
#' library(VariantAnnotation)
#'
#' # load default SEMs
#' data(SEMC)
#'
#' # create a VRanges object
#' vr <- VRanges(
#'     seqnames = "chr12",
#'     ranges = 94136009,
#'     ref = "G", alt = "C"
#' )
#'
#' # calculate binding propensity
#' s <- scoreVariants(vr, SEMC, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
#'
#' getRanges(s)
#'
setMethod(
    "getRanges", "SEMScores",
    function(x) x@ranges
)


#' Accessor semData slot in a SEMScores object
#'
#' @param x a SEMScores object
#' @rdname semData
#' @export
#'
#' @examples
#' library(VariantAnnotation)
#'
#' # load default SEMs
#' data(SEMC)
#'
#' # create a VRanges object
#' vr <- VRanges(
#'     seqnames = "chr12",
#'     ranges = 94136009,
#'     ref = "G", alt = "C"
#' )
#'
#' # calculate binding propensity
#' s <- scoreVariants(vr, SEMC, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
#'
#' semData(s)
#'
setMethod(
    "semData", "SEMScores",
    function(x) x@semData
)


#' Accessor scores slot in a SEMScores object
#'
#' @param x a SEMScores object
#'
#' @rdname scores
#' @keywords internal
#' @export
#'
#' @examples
#' library(VariantAnnotation)
#'
#' # load default SEMs
#' data(SEMC)
#'
#' # create a VRanges object
#' vr <- VRanges(
#'     seqnames = "chr12",
#'     ranges = 94136009,
#'     ref = "G", alt = "C"
#' )
#'
#' # calculate binding propensity
#' s <- scoreVariants(vr, SEMC, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
#'
#' scores(s)
#'
setMethod(
    "scores", "SEMScores",
    function(x) x@scores
)

setMethod("scores<-", "SEMScores", function(x, value) {
    x@scores <- value
    x
})


# ---- show ----


#' Show method for SEMScores objects
#'
#' Prints information about the number of variants, SEM meta data columns, and
#' the scoring table if scoreVariants has been run.
#'
#' @param object a SEMScores object
#'
#' @importFrom methods show
#'
#' @return An invisible NULL
#'
#' @rdname show-SEMScores
#'
#' @export
setMethod(
    "show", "SEMScores",
    function(object) {
        cat("An object of class SEMScores\n")

        # show ranges
        num_vars <- length(object@ranges)
        vars_id_list <- .formatList(x = object@ranges$id)

        cat("ranges(", num_vars, "): ", sep = "")
        cat(paste(vars_id_list, collapse = " "))

        # show semData
        meta_cols <- names(object@semData)
        n_meta_cols <- length(meta_cols)
        meta_cols_list <- .formatList(x = meta_cols)

        cat("\nsemData(", n_meta_cols, "): ", meta_cols_list, sep = "")

        # show scores
        n_scores <- nrow(object@scores)
        cat("\nscores(", n_scores, "):\n", sep = "")
        if (n_scores > 0) {
            print(object@scores)
        }
    }
)
