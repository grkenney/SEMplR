# ---- constructor ----

#' SNPEffectMatrix object and constructor
#'
#' Constructs a SNPEffectMatrix class object.
#'
#' @param sem A `data.table` object of format base position (rows) by nucleic
#' acid (columns). Column names must inclue A, C, T, and G; other columns
#' will be ignored.
#' @param baseline A `numeric` scrambled baseline, representing the binding
#' score of randomly scrambled kmers of the same length.
#' @param semId A `character` unique identifier for the SEM.
#' @importFrom methods new
#'
#' @return a SNPEffectMatrix object
#' @docType class
#' @export
#'
#' @examples
#' # build a matrix for a motif of length 4
#' m <- matrix(rnorm(16), nrow = 4)
#' colnames(m) <- c("A", "C", "G", "T")
#'
#' # build a SNPEffectMatrix object
#' SNPEffectMatrix(m, baseline = 1, semId = "sem_id")
#'
SNPEffectMatrix <- function(sem, baseline, semId) {
    # convert matrix to data.table in case given in another format
    sem <- data.table::as.data.table(sem)

    # drop columns not names after nucleotides
    sem <- sem[, c("A", "C", "G", "T")]

    # convert baseline to numeric
    baseline <- as.numeric(baseline)

    new("SNPEffectMatrix",
        sem = sem,
        baseline = baseline,
        semId = semId
    )
}


# ---- accessors ----

#' Access the SEM from a SNPEffectMatrix object
#'
#' @param x SNPEffectMatrix object
#'
#' @returns A position (rows) by nucleic acid (columns)
#' `data.table` is returned
#' @examples
#' # Isolate a single SNPEffectMatrix object from the default 
#' # SNPEffectMatrixCollection
#' sm <- sems(sc)[[1]]
#' 
#' # Access the matrix
#' getSEM(sm)
#'
#' @rdname getSEM
#' @export
setMethod(
    "getSEM", "SNPEffectMatrix",
    function(x) x@sem
)


#' Access baseline from a SNPEffectMatrix object
#'
#' @param x SNPEffectMatrix object
#'
#' @returns The `numeric` baseline value is returned
#'
#' @examples
#' # Isolate a single SNPEffectMatrix object from the default 
#' # SNPEffectMatrixCollection
#' sm <- sems(sc)[[1]]
#' 
#' # Access the baseline
#' getBaseline(sm)
#' 
#' @rdname getBaseline
#' @export
setMethod(
    "getBaseline", "SNPEffectMatrix",
    function(x) x@baseline
)


#' Access semId from a SNPEffectMatrix object
#'
#' @param x SNPEffectMatrix object
#'
#' @returns The `character` id is returned
#'
#' @examples
#' # Isolate a single SNPEffectMatrix object from the default 
#' # SNPEffectMatrixCollection
#' sm <- sems(sc)[[1]]
#' 
#' # Access the SEM id
#' getSEMId(sm)
#'
#' @rdname getSEMId
#' @export
setMethod(
    "getSEMId", "SNPEffectMatrix",
    function(x) x@semId
)


# ---- show ----

#' show for SNPEffectMatrix
#'
#' @param object SNPEffectMatrix object.
#'
#' @importFrom methods show 
#'
#' @rdname SNPEffectMatrix-class
setMethod(
    "show", "SNPEffectMatrix",
    function(object) {
        cat("An object of class SNPEffectMatrix\n")
        cat("semId: ", object@semId)
        cat("\nbaseline: ", object@baseline)
        cat("\nsem:\n")
        print(object@sem)
    }
)
