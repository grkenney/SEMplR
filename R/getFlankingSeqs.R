#' Get up and downstream sequences of a variant
#'
#' Finds nucleotides up and downstream of a variant within a genome
#'
#' @param x A `VRanges` or `GRanges` object with one or more variants. 
#' `seqnames` and `ranges` fields are required.
#' @param up Numeric, number of bases to return upstream of variant
#' @param down Numeric, number of bases to return downstream of variant
#' @param bs_genome_obj A `BSgenome` object for the genome build to use.
#' Defaults to `BSgenome.Hsapiens.UCSC.hg19`.
#' @param variant Boolean to indicate whether the provided query is a variant
#' or a single position.
#' @param allele If providing a `GRanges` object, can optionally provide an allele
#' to score for the indicated range. Otherwise, will query the BSgenome object
#' for the sequence at the indicated genomic position in x.
#'
#' @return a named nested list with`upstream`, `downstream`,
#' `ref_seq`, and `alt_seq` fields when `x` is a `VRanges` object or just 
#' `upstream`, `downstream`, and `seq` if `x` is a `GRanges` object.
#'
querySeqs <- function(x, up, down,
                      bs_genome_obj,
                      variant=T,
                      allele=NULL){
  # get start position of variant
  start_pos <- BiocGenerics::start(x)
  end_pos <- BiocGenerics::end(x)

  # check that ref allele matches BSgenome ref
  bsgenome_ref <- Biostrings::getSeq(bs_genome_obj,
                                     names=GenomeInfoDb::seqnames(x),
                                     start=start_pos,
                                     end=start_pos) |> as.character()

  flanking_seqs <- Biostrings::getSeq(bs_genome_obj,
                                      names=GenomeInfoDb::seqnames(x),
                                      start=start_pos-up,
                                      end=end_pos+down)
  
  if (variant) {
    # if the BSgenome reference allele does not match the provided reference
    # allele, give a warning
    # if (bsgenome_ref != VariantAnnotation::ref(x)) {
    #   warning(sprintf("Provided ref '%s' != BSgenome ref '%s' for variant at %s:%d \n",
    #                   VariantAnnotation::ref(x),
    #                   bsgenome_ref,
    #                   GenomeInfoDb::seqnames(x),
    #                   start_pos))
    # }
    
    concat_ref_seq <- Biostrings::xscat(flanking_seqs[1:up],
                                        as.character(VariantAnnotation::ref(x)),
                                        flanking_seqs[(up+2):(up+1+down)])
    concat_alt_seq <- Biostrings::xscat(flanking_seqs[1:up],
                                        as.character(VariantAnnotation::alt(x)),
                                        flanking_seqs[(up+2):(up+1+down)])
    # combine all data into a named vector
    seq_data <- c(upstream = as.character(flanking_seqs[1:up]),
                  downstream = as.character(flanking_seqs[(up+2):(up+1+down)]),
                  ref_seq = as.character(concat_ref_seq),
                  alt_seq = as.character(concat_alt_seq))
  } else {
    if (is.null(allele)) {
      allele <- bsgenome_ref
    }
    concat_seq <- Biostrings::xscat(flanking_seqs[1:up],
                                    allele,
                                    flanking_seqs[(up+2):(up+1+down)])
    # combine all data into a named vector
    seq_data <- c(upstream = as.character(flanking_seqs[1:up]),
                  downstream = as.character(flanking_seqs[(up+2):(up+1+down)]),
                  seq = as.character(concat_seq))
  }
  return(seq_data)
}


#' Get up and downstream sequences of variant(s) and store in VRanges metadata
#'
#' Finds nucleotides up and downstream of a variant within a genome and store
#' within `upstream` and `downstream` fields in the metadata of the supplied
#' `VRanges` object.
#'
#' @param x A `VRanges` or `GRanges` object with one or more variants. 
#' seqnames and ranges fields are required.
#' @param up Numeric, number of bases to return upstream of variant
#' @param down Numeric, number of bases to return downstream of variant
#' @param bs_genome_obj A `BSgenome` object for the genome build to use.
#' Defaults to `BSgenome.Hsapiens.UCSC.hg19`.
#'
#' @importFrom methods is
#'
#' @return a `VRanges` object with metadata columns `upstream` and `downstream`.
#'
#' @export
getFlankingSeqs <- function(x, up, down,
                            bs_genome_obj) {
  if (is(x)[1] != "VRanges" & is(x)[1] != "GRanges") {
    stop("Input argument x must be of class VRanges or GRanges")
  }
  
  if (length(x) < 1) {
    stop("VRanges object must contain at least one variant.")
  }
  
  # define metadata columns in x to store flanking seqs
  if(is(x)[1] == "VRanges") {
    meta_cols <- c("upstream", "downstream", "ref_seq", "alt_seq")
    variant <- T
  } else {
    meta_cols <- c("upstream", "downstream", "seq")
    variant <- F
  }
  S4Vectors::mcols(x)[meta_cols] <- NA

  # iterate over variants (rows) in vrange object
  flanking_seqs <- lapply(seq_along(x),
                          function(i) querySeqs(x[i],
                                                up = up, down = down,
                                                bs_genome_obj = bs_genome_obj, 
                                                variant = variant))

  # add flanking sequence to vranges metadata
  S4Vectors::mcols(x)[meta_cols] <- as.data.frame(do.call(rbind, 
                                                          flanking_seqs))
  return(x)
}
