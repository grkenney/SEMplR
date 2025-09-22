# ---- constructor ----

#' SNPEffectMatrixCollection object and constructor
#'
#' Constructs a SNPEffectMatrixCollection class object.
#'
#' @param sems A list of SNPEffectMatrix objects
#' @param semData A `data.table` containing meta data for each SEM
#' @param semKey A column name in semData corresponding to the semIds in each 
#' SNPEffectMatrix object in the sems parameter
#' 
#' @importFrom methods new
#'
#' @return a SNPEffectMatrixCollection object
#' @docType class
#' @export
#' 
#' @examples
#' # 
#' # build two matries for a motif of length 4
#' m1 <- matrix(rnorm(16), nrow = 4)
#' m2 <- matrix(rnorm(16), nrow = 4)
#' colnames(m1) <- c("A", "C", "G", "T")
#' colnames(m2) <- c("A", "C", "G", "T")
#' 
#' # build two SNPEffectMatrix objects
#' s1 <- SNPEffectMatrix(m1, baseline = 1, semId = "sem_id_1")
#' s2 <- SNPEffectMatrix(m2, baseline = 1, semId = "sem_id_2")
#' 
#' # build a meta data table
#' md <- data.table::data.table(transcription_factor = c("tf1", "tf2"), 
#'                              cell_type = c("HepG2", "K562"), 
#'                              sem_id = c("sem_id_1", "sem_id_2"))
#' 
#' # build a SNPEffectMatrixCollection object
#' SNPEffectMatrixCollection(list(s1, s2), semData = md, semKey = "sem_id")
#' 
SNPEffectMatrixCollection <- function(sems, semData = NULL, semKey = "") {
  SEM_KEY <- .SD <- NULL
  
  # if a single SNPEffectMatrix object, make it a list
  if (is(sems, "SNPEffectMatrix")) {
    sems <- list(sems)
  }
  
  # use semIds for list names
  names(sems) <- lapply(sems, function(x) unlist(x@semId))
  
  # must supply a key if supplying semData
  if (!is.null(semData) & semKey == "") {
    stop("must provide a semKey if providing semData",
    " See ?SNPEffectMatrixCollection")
  }
  
  # convert semData to data.table
  semData <- data.table::data.table(semData)
  
  if (nrow(semData) > 0) {
    # check that semKey is a column in semData
    if (!(semKey %in% names(semData))) {
      stop("semKey must be a column in semData \n\t '", semKey, 
           "' is not in names(semData)")
    }
    
    # if semKey column has a ".sem" suffix, remove
    # set the data.table key to the semKey column
    if (any(grepl(".sem", semData[, .SD, .SDcols = semKey]))) {
      semData[, SEM_KEY := lapply(semData[, .SD, .SDcols = semKey], 
                                  \(x) gsub(".sem", "", x))]
      data.table::setkey(semData, SEM_KEY)
      semKey <- "SEM_KEY"    # update key
      rlang::inform("Removing .sem suffixes from semKey.")
      message("formatted key now stored in column 'SEM_KEY'.")
    } else {
      data.table::setkeyv(semData, semKey)
    }
  }
  
  new("SNPEffectMatrixCollection",
      sems = sems,
      semData = semData,
      semKey = semKey
  )
}


# ---- accessors ----

#' Access sems matrix from a SNPEffectMatrixCollection object
#'
#' @param x SNPEffectMatrixCollection object
#' @param semId optional `character` corresponding to an SEM in the 
#' SNPEffectMatrixCollection object. See semIds with `names(sems(x))`. 
#' Defaults to returning all SEMs.
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
#' sc <- SNPEffectMatrixCollection(s)
#'
#' ## Access count matrix
#' sems(sc)
#'
#' @rdname sems
#' @export
setMethod("sems", "SNPEffectMatrixCollection", 
          function(x, semId = NULL) {
            if (is.null(semId)) {
              x@sems
            } else if (length(semId) == 1) {
              x@sems[[semId]]
            } else {
              x@sems[semId]
            }
          })


#' Access semData from a SNPEffectMatrixCollection object
#'
#' @param x SNPEffectMatrixCollection object
#'
#' @returns A `data.table` is returned
#'
#' @examples
#' ## Create example SEM
#' df <- data.frame(A=c(1, 2, 3),
#'                  C=c(1, 2, 3),
#'                  G=c(1, 2, 3),
#'                  T=c(1, 2, 3))
#' sem_data <- data.frame(SEM = "motif_id", TF = "tf_id")
#' s <- SNPEffectMatrix(df, 1.205, "motif_id")
#' sc <- SNPEffectMatrixCollection(s, semData = sem_data, semKey = "SEM")
#' 
#' ## Access count matrix
#' semData(sc)
#'
#' @rdname semData
#' @export
setMethod("semData", "SNPEffectMatrixCollection", 
          function(x) x@semData)


# ---- show ----

#' Show method for SNPEffectMatrixCollection objects
#'
#' Prints information about the number of SEMs and meta data columns included
#' in the object.
#'
#' @param object a SNPEffectMatrixCollection object
#' 
#' @importFrom methods show
#' 
#' @rdname show
setMethod("show", "SNPEffectMatrixCollection",
          function(object) {
            cat("An object of class SNPEffectMatrixCollection\n")
            num_sems <- length(object@sems)
            
            s <- lapply(object@sems, function(x) x@semId) |>
              unlist()
            sem_id_list <- .formatList(s)
            
            cat("sems(",num_sems,"): ", sep = "")
            cat(paste(sem_id_list, collapse = " "))
            
            meta_cols <- names(object@semData)
            n_cols <- length(meta_cols)
            meta_cols_list <- .formatList(meta_cols)
            
            cat("\nsemData(", n_cols, "): ", meta_cols_list, sep = "")
          }
)

