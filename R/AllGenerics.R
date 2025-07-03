# ---- SNPEffectMatrix ----

#' @rdname sem
#' @export
setGeneric("sem", 
           function(x) standardGeneric("sem"))

#' @rdname baseline
#' @export
setGeneric("baseline", 
           function(x) standardGeneric("baseline"))

#' @rdname semId
#' @export
setGeneric("semId", 
           function(x) standardGeneric("semId"))


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

#' @rdname variants
#' @export
setGeneric("variants", function(x) standardGeneric("variants"))

#' @rdname scores
#' @export
setGeneric("scores", function(x) standardGeneric("scores"))
setGeneric("scores<-", 
           function(x, value) standardGeneric("scores<-"))
