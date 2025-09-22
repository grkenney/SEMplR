# format a list of items to print well in show function
.formatList <- function(x) {
  num_items <- length(x)
  if (num_items > 5) {
    first2 <- x[seq_len(2)]
    last2 <- x[(num_items-1):num_items]
    formatted_list <- paste0(paste(first2, collapse = ", "), " ... ", 
                             paste(last2, collapse = ", "))
  } else {
    formatted_list <- x |>
      paste0(collapse = ", ")
  }
  return(formatted_list)
}


# most functions in this package expect a SNPEffectMatrixCollection
# for convenience, this functions converts lists of SNPEffectMatrix's,
# or single SNPEffectMatrix into SNPEffectMatrixCollections
.convertToSNPEffectMatrixCollection <- function(x) {
  # if it's already a collection, return as is
  # if given a list of SNPEffectMatrices or a single SNPEffectMatrix, 
  # make a collection
  # else, fail if unable to convert
  if (is(x, "SNPEffectMatrixCollection")) {
    return(x)
  } else if(is(x, "list") | is(x, "vector")) {
    # check that all elements are SNPEffectMatrices
    class_check <- lapply(x, function(y) 
      is(y, class2 = "SNPEffectMatrix")) |>
      unlist()
    if (all(class_check)) {
      return(SNPEffectMatrixCollection(x))
    } else {
      invalid_class <- is(x[!class_check][1])[1]
      rlang::abort(paste0(
        "unable to convert object of class ", invalid_class, 
        " to class SNPEffectMatrixCollection.\n",
        "See ?SNPEffectMatrixCollection or use the provided default 'sc'"))
    }
    
  } else if (is(x, "SNPEffectMatrix")) {
    return(SNPEffectMatrixCollection(sems = x))
  } else {
    rlang::abort(paste0(
      "unable to convert object of class ", is(x)[1], 
      " to class SNPEffectMatrixCollection.\n",
      "See ?SNPEffectMatrixCollection or use the provided default 'sc'"))
  }
}


# Given a single VRange or GRange, construct a unique id from the 
# position and allele information
.makeVariantId <- \(x, refCol = NULL, altCol = NULL) {
  start_pos <- IRanges::start(IRanges::ranges(x))
  end_pos <- IRanges::end(IRanges::ranges(x))
  sn <- GenomeInfoDb::seqnames(x)
  
  if (is(x, "VRanges")) {
    ref_allele <- as.character(VariantAnnotation::ref(x))
    alt_allele <- as.character(VariantAnnotation::alt(x))
  } else {
    if (is.null(refCol) | is.null(altCol)) {
      rlang::abort(paste0("If providing a GRanges object, ",
                          "both refCol and altCol must be defined"))
    }
    ref_allele <- as.character(S4Vectors::mcols(x)[, refCol])
    alt_allele <- as.character(S4Vectors::mcols(x)[, altCol])
  }
  
  if (ref_allele == "") {
    allele_str <- paste0("ins", alt_allele)
  } else if (alt_allele == "") {
    allele_str <- paste0("del", ref_allele)
  } else {
    allele_str <- paste0(ref_allele, ">", alt_allele)
  }
  
  pos_str <- ifelse(start_pos == end_pos, 
                    start_pos, paste0(start_pos, "-", end_pos))
  
  vid <- paste0(sn, ":", pos_str, ":", allele_str)
  return(vid)
}


# convert a SNP Effect Matrix to a Position Probability Matrix
.semToPpm <- \(s) {
  # normalize matrix
  norm_score <- apply(getSEM(s), 1, 
                      function(x) (2^x - 2^getBaseline(s)) / 
                        abs(2^getBaseline(s)))
  # replace negative scores with zero
  norm_score[norm_score < 0] <- 0
  # make all rows sum to 1
  ppm <- apply(norm_score, 2, function(x) x / sum(x))
  return(ppm)
}
