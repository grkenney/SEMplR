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
SNPEffectMatrix <- function(sem, baseline, semId) {
  .SD <- NULL
  
  # convert matrix to data.table in case given in another format
  sem <- data.table::as.data.table(sem)
  
  # drop first column
  sem <- sem[, .SD, .SDcols = c("A", "C", "G", "T")]
  
  # convert baseline to numeric
  baseline <- as.numeric(baseline)
  
  new("SNPEffectMatrix",
      sem = sem,
      baseline = baseline,
      semId = semId
  )
}


# ---- accessors ----

#' Access sem matrix from a SNPEffectMatrix object
#'
#' @param x SNPEffectMatrix object
#'
#' @returns A position (rows) by nucleic acid (columns) 
#' `data.table` is returned
#'
#' @examples
#' ## Create example SEM
#' df <- data.frame(A=c(1, 2, 3),
#'                  C=c(1, 2, 3),
#'                  G=c(1, 2, 3),
#'                  T=c(1, 2, 3))
#' s <- SNPEffectMatrix(df, 1.205, "motif_id")
#'
#' ## Access count matrix
#' sem(s)
#'
#' @rdname sem
#' @export
setMethod("sem", "SNPEffectMatrix", 
          function(x) x@sem)


#' Access baseline from a SNPEffectMatrix object
#'
#' @param x SNPEffectMatrix object
#'
#' @returns The `numeric` baseline value is returned
#'
#' @examples
#' ## Create example SEM
#' df <- data.frame(A=c(1, 2, 3),
#'                  C=c(1, 2, 3),
#'                  G=c(1, 2, 3),
#'                  T=c(1, 2, 3))
#' s <- SNPEffectMatrix(df, 1.205, "motif_id")
#'
#' ## Access count matrix
#' baseline(s)
#'
#' @rdname baseline
#' @export
setMethod("baseline", "SNPEffectMatrix", 
          function(x) x@baseline)


#' Access semId from a SNPEffectMatrix object
#'
#' @param x SNPEffectMatrix object
#'
#' @returns The `character` semId is returned
#'
#' @examples
#' ## Create example SEM
#' df <- data.frame(A=c(1, 2, 3),
#'                  C=c(1, 2, 3),
#'                  G=c(1, 2, 3),
#'                  T=c(1, 2, 3))
#' s <- SNPEffectMatrix(df, 1.205, "motif_id")
#'
#' ## Access count matrix
#' semId(s)
#'
#' @rdname semId
#' @export
setMethod("semId", "SNPEffectMatrix", 
          function(x) x@semId)


# ---- show ----

#' show for SNPEffectMatrix
#' 
#' @param object SNPEffectMatrix object.
#' @rdname SNPEffectMatrix-class
#' 
#' @export
setMethod("show", "SNPEffectMatrix",
          function(object) {
            cat("An object of class SNPEffectMatrix\n")
            cat("semId: ", object@semId)
            cat("\nbaseline: ", object@baseline)
            cat("\nsem:\n")
            print(object@sem)
          }
)