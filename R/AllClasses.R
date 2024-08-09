
## SemplR class ----------------------------------------------------------------

#' Class for storing SEM motif binding calculations for multiple variants
#'
#' @slot variants A `VRanges` object to hold one or more variants
#' @slot scores A `data.table` object for motif information and binding scores
#'
#' @export
setClass("SemplR",
         slots = c(variants = "VRanges",
                   scores = "data.table"
                   ),
         prototype = list(
           variants = VRanges(),
           scores = data.frame()
         )
)

#' SemplR object and constructor
#'
#' Constructs a SemplR class object.
#'
#' @param variants A `VRanges` object to hold one or more variants
#' @param scores A `data.table` object for motif information and binding scores
#'
#' @return a SemplR object
#' @docType class
#' @aliases SemplR-class
#' @rdname SemplR
#' @export
SemplR <- function(variants, scores) {
  new("SemplR",
      variants = VRanges(),
      scores = data.table()
      )
}
