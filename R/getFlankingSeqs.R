#' Get up and downstream sequences of a variant
#'
#' Finds nucleotides up and downstream of a variant within a genome
#'
#' @param vr A `VRanges` object with one or more variants. seqnames and ranges
#' fields are required.
#' @param up Numeric, number of bases to return upstream of variant
#' @param down Numeric, number of bases to return downstream of variant
#' @param bs_genome_obj A `BSgenome` object for the genome build to use.
#' Defaults to `BSgenome.Hsapiens.UCSC.hg19`.
#'
#' @return a named nested list with`upstream` and `downstream`,
#' `ref_seq`, and `alt_seq` fields.
#'
query_seqs <- function(vr, up, down,
                       bs_genome_obj=BSgenome.Hsapiens.UCSC.hg19::Hsapiens){
  # get start position of variant
  start_pos <- BiocGenerics::start(vr)
  end_pos <- BiocGenerics::end(vr)

  # check that ref allele matches BSgenome ref
  bsgenome_ref <- Biostrings::getSeq(bs_genome_obj,
                                     names=GenomeInfoDb::seqnames(vr),
                                     start=start_pos,
                                     end=start_pos) |> as.character()

  # if the BSgenome reference allele does not match the provided reference
  # allele, give a warning
  # if (bsgenome_ref != VariantAnnotation::ref(vr)) {
  #   warning(sprintf("Provided ref '%s' != BSgenome ref '%s' for variant at %s:%d \n",
  #                   VariantAnnotation::ref(vr),
  #                   bsgenome_ref,
  #                   GenomeInfoDb::seqnames(vr),
  #                   start_pos))
  # }

  flanking_seqs <- Biostrings::getSeq(bs_genome_obj,
                                      names=GenomeInfoDb::seqnames(vr),
                                      start=start_pos-up,
                                      end=end_pos+down)
  
  concat_ref_seq <- Biostrings::xscat(flanking_seqs[1:up],
                    as.character(VariantAnnotation::ref(vr)),
                    flanking_seqs[(up+2):(up+1+down)])
  concat_alt_seq <- Biostrings::xscat(flanking_seqs[1:up],
                                      as.character(VariantAnnotation::alt(vr)),
                                      flanking_seqs[(up+2):(up+1+down)])

  # combine all data into a named vector
  seq_data <- c(upstream = as.character(flanking_seqs[1:up]),
                downstream = as.character(flanking_seqs[(up+2):(up+1+down)]),
                ref_seq = as.character(concat_ref_seq),
                alt_seq = as.character(concat_alt_seq))
  return(seq_data)
}


#' Get up and downstream sequences of variant(s) and store in VRanges metadata
#'
#' Finds nucleotides up and downstream of a variant within a genome and store
#' within `upstream` and `downstream` fields in the metadata of the supplied
#' `VRanges` object.
#'
#' @param vr A `VRanges` object with one or more variants. seqnames and ranges
#' fields are required.
#' @param up Numeric, number of bases to return upstream of variant
#' @param down Numeric, number of bases to return downstream of variant
#' @param bs_genome_obj A `BSgenome` object for the genome build to use.
#' Defaults to `BSgenome.Hsapiens.UCSC.hg19`.
#'
#' @return a `VRanges` object with metadata columns `upstream` and `downstream`.
#'
#' @export
getFlankingSeqs <- function(vr, up, down,
                            bs_genome_obj=BSgenome.Hsapiens.UCSC.hg19::Hsapiens) {
  if (length(vr) < 1) {
    stop("Vranges object must contain at least one variant.")
  }
  
  # initialize metadata columns in vrange to store up and downstream seqs
  vr$upstream <- NA
  vr$downstream <- NA

  vr$ref_seq <- NA
  vr$alt_seq <- NA

  # iterate over variants (rows) in vrange object
  flanking_seqs <- lapply(seq_along(vr),
                          function(i) query_seqs(vr[i],
                                                 up = up, down = down,
                                                 bs_genome_obj = bs_genome_obj))
  # add flanking sequence to vranges metadata
  S4Vectors::mcols(vr)[c("upstream", "downstream", 
              "ref_seq", "alt_seq")] <- as.data.frame(do.call(rbind, 
                                                              flanking_seqs))

  return(vr)
}
