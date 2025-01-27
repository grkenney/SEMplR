#' View the top scoring frames for a given variant and SEM
#'
#' @param sempl_obj SemplScores object
#' @param score_index index of score to view
#'
#' @return a SequenceFrame object
#'
#' @export
viewFrames <- function(sempl_obj, score_index) {
  v <- variants(sempl_obj)[variants(sempl_obj)$id == scores(sempl_obj)$varId[score_index]]
  
  if (ref(v) == "") {
    vref <- rep("-", nchar(alt(v)))
  } else {
    vref <- ref(v)
  }
  
  if (alt(v) == "") {
    valt <- rep("-", nchar(ref(v)))
  } else {
    valt <- alt(v)
  }
  
  rs <- paste0(v$upstream,
               vref,
               v$downstream)
  vs <- paste0(v$upstream,
               valt,
               v$downstream)
  
  vi <- (nchar(v$upstream)+1):(nchar(v$upstream) + nchar(vref))

  sf <- SequenceFrame(c("ref", "alt"), 
                      sequence = c(rs, vs),
                      frameStart = c(sempl_obj@scores$nonRiskVarIndex[score_index],
                                     sempl_obj@scores$riskVarIndex[score_index]), 
                      frameStop = c(sempl_obj@scores$nonRiskVarIndex[score_index] +
                                      nrow(sems[[sempl_obj@scores$semId[score_index]]]@sem),
                                    sempl_obj@scores$riskVarIndex[score_index] +
                                      nrow(sems[[sempl_obj@scores$semId[score_index]]]@sem)),
                      variantIndex = vi)
  return(sf)
}
