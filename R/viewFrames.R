#' View the top scoring frames for a given variant and SEM
#'
#' @param s SemplScores object
#' @param vid variant id corresponding to the frame to view
#' @param sid SEM id corresponding to the frame to view
#' @param scoreIndex index of score to view. `vid` and `sid` takes prescedent 
#' over this parameter.
#'
#' @return a SequenceFrame object
#'
#' @export
viewFrames <- function(s, vid = NULL, sid = NULL, scoreIndex = NULL) {
  varId <- semId <- NULL
  
  # get the corresponding sempl scores row given the vid and sid or the index 
  if (!is.null(vid) & !is.null(sid)) {
    ss <- scores(s)[varId == vid & semId == sid]
    
    # stop if no rows found matching vid and sid
    if (nrow(ss) == 0) {
      stop("no entry matching vid == ", vid, " and sid == ", sid, " found in",
           " the score table.")
    }
    
  } else if (!is.null(scoreIndex)) {
    # check that scoreIndex is within the range of s scores
    if (scoreIndex > nrow(scores(s))) {
      stop("scoreIndex ", scoreIndex, " out of range in SEMplScores object")
    }
    # get the corresponding sempl scores row
    ss <- scores(s)[scoreIndex]
    
  } else {
    # require that the user supply either a vid and sid or a scoreIndex
    stop("must provide either a scoreIndex or a vid and sid parameter",
         " corresponding to the frame to view.")
  }
  
  # get the variant that matches the score index
  vi <- getRanges(s)$id == ss$varId
  v <- getRanges(s)[vi]
  
  # get the corresponding sem matrix
  motifbp <- nchar(ss$refSeq)

  # get the reference and alternative alleles
  ref_allele <- as.character(VariantAnnotation::ref(v))
  alt_allele <- as.character(VariantAnnotation::alt(v))
  
  # if the ref allele is blank (deletion), repeat "-" for length of alt allele
  if (ref_allele == "") {
    ref_str <- rep("-", nchar(alt_allele)) |> paste0(collapse = "")
  } else {
    ref_str <- ref_allele
  }
  
  # if the alt allele is blank (deletion), repeat "-" for length of ref allele
  if (alt_allele == "") {
    alt_str <- rep("-", nchar(ref_allele)) |> paste0(collapse = "")
  } else {
    alt_str <- alt_allele
  }
  
  # concatenate flanking sequences to allele strings
  full_ref_str <- paste0(v$upstream,
                         ref_str,
                         v$downstream)
  full_alt_str <- paste0(v$upstream,
                         alt_str,
                         v$downstream)
  
  # get the index of the variant location on the sequence
  vi <- (nchar(v$upstream)+1):(nchar(v$upstream) + nchar(ref_str))

  # get the frame start and stop index
  fstart <- c(ss$refVarIndex,
              ss$altVarIndex)
  fstop <- c(fstart[1] + motifbp-1,
             fstart[2] + motifbp-1)
  
  # check if fstart is before the first basepair or if fstop is after the 
  # last basepair
  if (any(fstart < 1) | any(fstop > nchar(full_ref_str))) {
    stop("Frame index range (", min(fstart), ":", max(fstop), ") exceeds ",
         "sequence length (", nchar(full_ref_str),")")
  }
  
  # warn if the variant index is outside the frame index range for ref and alt
  if (!any(vi %in% fstart[1]:fstop[1]) | !any(vi %in% fstart[2]:fstop[2])) {
    warning("Variant index ", paste0(vi, collapse = ", "), " is not within ref",
    "sequence frame.")
  }
  
  # create the SequenceFrame object
  sf <- SequenceFrame(c("ref", "alt"), 
                      sequence = c(full_ref_str, full_alt_str),
                      frameStart = fstart, 
                      frameStop = fstop,
                      variantIndex = vi,
                      motif = ss$semId,
                      variantName = as.character(ss$varId))
  return(sf)
}
