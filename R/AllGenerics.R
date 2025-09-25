# ---- SNPEffectMatrix ----

#' @rdname getSEM
#' @export
setGeneric(
    "getSEM",
    function(x) standardGeneric("getSEM")
)

#' @rdname getBaseline
#' @export
setGeneric(
    "getBaseline",
    function(x) standardGeneric("getBaseline")
)

#' @rdname getSEMId
#' @export
setGeneric(
    "getSEMId",
    function(x) standardGeneric("getSEMId")
)


# ---- SNPEffectMatrixCollection ----

#' @rdname getSEMs
#' @export
setGeneric(
    "getSEMs",
    function(x, semId = NULL) standardGeneric("getSEMs")
)

#' @rdname semData
#' @export
setGeneric(
    "semData",
    function(x) standardGeneric("semData")
)


# ---- SEMplScores ----

#' @rdname getRanges
#' @export
setGeneric("getRanges", function(x) standardGeneric("getRanges"))

#' @rdname scores
#' @export
setGeneric("scores", function(x) standardGeneric("scores"))
setGeneric(
    "scores<-",
    function(x, value) standardGeneric("scores<-")
)
