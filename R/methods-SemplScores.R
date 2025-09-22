# ---- constructor ----

#' SEMplScores object and constructor
#'
#' Constructs a SEMplScores class object.
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
#' @return a SEMplScores object
#' @docType class
#' @export
#' 
#' @examples
#' # load default SEMs
#' data(sc)
#' 
#' # create a VRanges object
#' vr <- VariantAnnotation::VRanges(seqnames = c("chr12", "chr19"),
#'                                  ranges = c(94136009, 10640062), 
#'                                  ref = c("G", "T"), alt = c("C", "A"))
#' 
#' SEMplScores(ranges = vr, semData = semData(sc))
#' 
SEMplScores <- function(ranges=NULL, semData=NULL, scores=NULL) {
  # if no ranges provided, make an empty VRanges object
  if (all(is.null(ranges))) {
    r <-  VariantAnnotation::VRanges()
  } else {
    r <- ranges
  }
  
  # # if more than one VRange, and there is no id column, make a unique id from
  # # position information
  # if (length(vr) == 0) {
  #   S4Vectors::mcols(vr) <- data.frame(id=NA)
  # } else if (!("id" %in% names(S4Vectors::mcols(vr)))) {
  #   vid <- lapply(1:length(vr), \(i) .makeVariantId(vr = vr[i])) |>
  #     unlist()
  #   S4Vectors::mcols(vr)$id <- vid
  # }
  
  if (is.null(scores)){
    scores_table <- data.table()
  } else {
    scores_table <- scores
  }
  
  if (is.null(semData)){
    semData <- data.table()
  }
  
  new("SEMplScores",
      ranges = r,
      semData = semData,
      scores = scores_table
  )
}


# ---- accessors ----

#' Access ranges slot in a SEMplScores object
#' 
#' @param x a SEMplScores object
#' @rdname getRanges
#' @export
#' 
#' @examples
#' library(VariantAnnotation)
#' 
#' # load default SEMs
#' data(sc)
#' 
#' # create a VRanges object
#' vr <- VRanges(seqnames = "chr12",
#'               ranges = 94136009, 
#'               ref = "G", alt = "C")
#' 
#' # calculate binding propensity
#' s <- scoreVariants(vr, sc, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
#' 
#' getRanges(s)
#' 
setMethod("getRanges", "SEMplScores", 
          function(x) x@ranges)


#' Accessor semData slot in a SEMplScores object
#' 
#' @param x a SEMplScores object
#' @rdname semData
#' @export
#' 
#' @examples
#' library(VariantAnnotation)
#' 
#' # load default SEMs
#' data(sc)
#' 
#' # create a VRanges object
#' vr <- VRanges(seqnames = "chr12",
#'               ranges = 94136009, 
#'               ref = "G", alt = "C")
#' 
#' # calculate binding propensity
#' s <- scoreVariants(vr, sc, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
#' 
#' semData(s)
#' 
setMethod("semData", "SEMplScores", 
          function(x) x@semData)


#' Accessor scores slot in a SEMplScores object
#' 
#' @param x a SEMplScores object
#' 
#' @rdname scores
#' @keywords internal
#' @export
#' 
#' @examples
#' library(VariantAnnotation)
#' 
#' # load default SEMs
#' data(sc)
#' 
#' # create a VRanges object
#' vr <- VRanges(seqnames = "chr12",
#'               ranges = 94136009, 
#'               ref = "G", alt = "C")
#' 
#' # calculate binding propensity
#' s <- scoreVariants(vr, sc, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
#' 
#' scores(s)
#' 
setMethod("scores", "SEMplScores", 
          function(x) x@scores)

setMethod("scores<-", "SEMplScores", function(x, value) {
  x@scores <- value
  x
})


# ---- show ----


#' Show method for SEMplScores objects
#'
#' Prints information about the number of variants, SEM meta data columns, and
#' the scoring table if scoreVariants has been run.
#'
#' @param object a SEMplScores object
#' 
#' @importFrom methods show
#' 
#' @rdname show
#' 
#' @export
setMethod("show", "SEMplScores",
          function(object) {
            cat("An object of class SEMplScores\n")
            
            # show ranges
            num_vars <- length(object@ranges)
            vars_id_list <- .formatList(x = object@ranges$id)
            
            cat("ranges(",num_vars,"): ", sep = "")
            cat(paste(vars_id_list, collapse = " "))
            
            # show semData
            meta_cols <- names(object@semData)
            n_meta_cols <- length(meta_cols)
            meta_cols_list <- .formatList(x = meta_cols)

            cat("\nsemData(", n_meta_cols, "): ", meta_cols_list, sep = "")
            
            # show scores
            n_scores <- nrow(object@scores)
            cat("\nscores(", n_scores, "):\n", sep="")
            if(n_scores > 0) { print(object@scores) }
          }
)


