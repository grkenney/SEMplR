# query sequence associates with a given range, x.
# b: BSgenome object, sp: starting position, ep: ending position
# sp and ep will override the start and end sites defined in range x
.querySeq <- function(x, b, sp = NULL, ep = NULL) {
  # get position info
  chromosome <- GenomeInfoDb::seqnames(x)

  if (is.null(sp)) {
    sp <- BiocGenerics::start(x)
  }
  
  if (is.null(ep)) {
    ep <- BiocGenerics::end(x)
  }
  
  # get sequence for given range
  qs <- Biostrings::getSeq(b, 
                           names = GenomeInfoDb::seqnames(x), 
                           start = sp, end = ep) |>
    as.character()
  return(qs)
}


# check that allele parameter is within the querried sequence and that there
# are no more than 2 alleles given
# x: range being tested, qs: sequence from range in reference genome,
# allele: sequence provided, 
.checkAlleles <- function(x, qs, allele) {
  # check if qs is one of the alleles, if not, throw a warning
  if (!(qs %in% allele)) {
    rlang::warn(paste0("Allele '", paste(allele, collapse = "' or '"), 
                       "' does not match reference sequence '", qs, "' in ", x))
  }
  
  # stop if more than 2 alleles provided
  if (length(allele) > 2) {
    rlang::abort(paste0(
      "allele parameter '", paste(allele, collapse = "', '"), 
      "' for ", x, " is invalid.\n",
      "The number of alleles must be 2 or less."))
  }
}


# build sequence associated with range, x, with upstream and downstream flanks
# up: number of bps upstream, down: number of bps downstream, b: BSgenome object
# allele: allele to replace the bps in the range
.buildSeq <- function(x, up, down, b, allele = NULL) {
    # query sequence of range given
    qs <- .querySeq(x, b)

    # query sequence of upstream flank
    sp_up <- BiocGenerics::start(x)-up
    ep_up <- BiocGenerics::start(x)-1

    qs_up <- .querySeq(x, b, sp = sp_up, ep = ep_up)
    
    # query sequence of downstream flank
    sp_down <- BiocGenerics::end(x)+1
    ep_down <- BiocGenerics::end(x)+down
    qs_down <- .querySeq(x, b, sp = sp_down, ep = ep_down)
    
    # if an allele is given, use it instead of qs
    if (is.null(allele)) {
      concat_seq <- paste0(qs_up, qs, qs_down)
      seq_data <- c(sequence = concat_seq)
    } else {
      .checkAlleles(x, qs, allele)
      
      concat_seq <- paste0(qs_up, allele, qs_down)
      seq_data <- c(ref_seq = concat_seq[1],
                    alt_seq = concat_seq[2])
    }
    
    return(seq_data)
}


# get allele information from a range object given an allele column name.
# if it's a VRange, use the built-in ref and alt alleles
.getAllele <- function(x, refCol = NULL, altCol = NULL) {
  # if a VRanges, use the built-in ref and alt alleles
  # else if GRanges, use the specified ref and alt col names to get alleles
  if(is(x, "VRanges")) {
    if (!is.null(refCol) | !is.null(altCol)) {
      rlang::inform(paste0(
        "VRanges object detected, ignoring allele column specifications.",
        " Alleles will be pulled using ref() and alt() functions."))
    }
    ra <- VariantAnnotation::ref(x) |> as.character()
    aa <- VariantAnnotation::alt(x) |> as.character()
    
  } else if(is(x, "GRanges")) {
    # stop if refCol or altCol are not valid column names
    invalid_cols <- c(refCol, altCol) %in% colnames(S4Vectors::mcols(x))
    if( !all(invalid_cols) ){
      invalid_col_names <- colnames(S4Vectors::mcols(x))[!invalid_cols]
      rlang::abort(paste0("'", invalid_col_names, 
                          "' is not a valid column name in mcols(x)"))
    }
    
    if (!is.null(refCol)) {
      ra <- S4Vectors::mcols(x)[refCol] |> unlist() |> unname()
    } else {
      ra <- NULL
    }
    
    if (!is.null(altCol)) {
      aa <- S4Vectors::mcols(x)[altCol] |> unlist() |> unname()
    } else {
      aa <- NULL
    }
  }
  
  alleles <- list(ref = ra, alt = aa)
  return(alleles)
}


#' Get sequence for genomic ranges and variants
#' 
#' Query sequences associated with a given range in a reference genome and 
#' optionally include nucleotides up and downstream of the range of interest
#' sequences. Store 
#' the results in the metadata of the range object.
#'
#' @param x A `VRanges` or `GRanges` object with one or more variants. 
#' seqnames and ranges fields are required.
#' @param genome A `BSgenome` object for the genome build to use.
#' @param up Numeric, number of bases to return upstream of variant
#' @param down Numeric, number of bases to return downstream of variant
#' @param refCol Column name in meta data storing the ref allele
#' @param altCol Column name in meta data storing the alt allele
#'
#' @importFrom methods is
#'
#' @return a `x` with added meta data columns. By default, only a "sequence"
#' column is added to the meta data. If variant information is supplied, either
#' in a `VRanges` object or through the "allele" parameter, "ref_seq" and 
#' "alt_seq" columns are added to the meta data.
#' 
#' @export
#' 
#' @examples
#' # set reference genome
#' b <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
#' 
#' # Query sequnece for a GRanges object
#' gr <- GenomicRanges::GRanges(seqnames = "chr12",
#'                              ranges = IRanges::IRanges(94136009, 94136019))
#' getRangeSeqs(gr, genome = b)
#' 
#' # Query variant sequences for a VRanges object with 10bp flanks
#' vr <- VariantAnnotation::VRanges(seqnames = "chr12",
#'                                  ranges = IRanges::IRanges(94136009),
#'                                  ref = "G", alt = "C")
#' getRangeSeqs(vr, genome = b, up = 10, down = 10)
#' 
getRangeSeqs <- function(x, genome, up = 0, down = 0,
                         refCol = NULL, altCol = NULL) {
  # stop if x is not a GRanges or VRanges object
  if ( !(is(x, "GRanges") | is(x, "VRanges")) ) {
    rlang::abort("x must be of class VRanges or GRanges")
  }
  
  # stop if x is empty
  if (length(x) < 1) {
    rlang::abort("x must contain at least one range")
  }
  
  a <- .getAllele(x = x, refCol = refCol, altCol = altCol)

  # iterate over ranges to find sequences
  flanking_seqs <- lapply(seq_along(x), 
                          function(i) .buildSeq(x[i], 
                                                up = up, 
                                                down = down, 
                                                b = genome,
                                                allele = c(a[["ref"]][i],
                                                           a[["alt"]][i])))
  flanking_df <- as.data.frame(do.call(rbind, flanking_seqs))
  # add flanking sequence to metadata
  meta_cols <- colnames(flanking_df)
  S4Vectors::mcols(x)[meta_cols] <- flanking_df
  return(x)
}
