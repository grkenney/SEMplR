# ---- show ----

#' Show function for SequenceFrame
#' 
#' @param object a SequenceFrame object
#' 
#' @importFrom crayon bgRed
#' @importFrom methods show
#' 
#' @rdname show
setMethod("show", "SequenceFrame",
          function(object) {
            for (i in 1:length(object@sequence)) {
              prefix <- substr(object@sequence[i], 
                               start = 1, 
                               stop = object@frameStart[i] - 1)
              
              fp <- substr(object@sequence[i], 
                           object@frameStart[i], 
                           min(object@variantIndex)-1)
              
              v <- substr(object@sequence[i], 
                          min(object@variantIndex),
                          max(object@variantIndex))
              
              fs <- substr(object@sequence[i], 
                           max(object@variantIndex)+1, 
                           object@frameStop[i])
              
              suffix <- substr(object@sequence[i], 
                               start = object@frameStop[i]+1, 
                               stop = nchar(object@sequence[i]))
              
              grey7 <- crayon::make_style("#777777", bg=TRUE)
              
              if (i == 1){
                cat("SEM:\t", object@motif, "\n")
                cat("Var:\t", object@variantName, "\n")
              }
              
              cat(object@seqName[i], ": ", prefix, grey7(fp), bgRed(v), grey7(fs), suffix, "\n", sep="")
            }
          }
)
