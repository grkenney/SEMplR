# ---- SNPEffectMatrix ----

#' @rdname getSEM
#' @export
setGeneric("getSEM", 
           function(x) standardGeneric("getSEM"))

#' @rdname getBaseline
#' @export
setGeneric("getBaseline", 
           function(x) standardGeneric("getBaseline"))

#' @rdname getSEMId
#' @export
setGeneric("getSEMId", 
           function(x) standardGeneric("getSEMId"))


# ---- SNPEffectMatrixCollection ----

#' @rdname sems
#' @export
setGeneric("sems", 
           function(x, semId = NULL) standardGeneric("sems"))

#' @rdname semData
#' @export
setGeneric("semData", 
           function(x) standardGeneric("semData"))


# ---- SemplScores ----

#' @rdname getRanges
#' @export
setGeneric("getRanges", function(x) standardGeneric("getRanges"))

#' @rdname scores
#' @export
setGeneric("scores", function(x) standardGeneric("scores"))
setGeneric("scores<-", 
           function(x, value) standardGeneric("scores<-"))
