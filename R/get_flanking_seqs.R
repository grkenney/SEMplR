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

  if (bsgenome_ref != VariantAnnotation::ref(vr)) {
    warning(sprintf("Provided ref '%s' != BSgenome ref '%s' for variant at %s:%d",
                    VariantAnnotation::ref(vr),
                    bsgenome_ref,
                    GenomeInfoDb::seqnames(vr),
                    start_pos))
  }

  # get sequence upstream and downstream of variant start/end
  upstream <- Biostrings::getSeq(bs_genome_obj,
                                 names=GenomeInfoDb::seqnames(vr),
                                 start=start_pos-up,
                                 end=start_pos-1)
  downstream <- Biostrings::getSeq(bs_genome_obj,
                                   names=GenomeInfoDb::seqnames(vr),
                                   start=end_pos+1,
                                   end=end_pos+down)

  # concatenate flanking sequences with variant and store in vrange metadata
  concat_ref_seq <- Biostrings::xscat(upstream,
                                      as.character(VariantAnnotation::ref(vr)),
                                      downstream)
  concat_alt_seq <- Biostrings::xscat(upstream,
                                      as.character(VariantAnnotation::alt(vr)),
                                      downstream)

  seq_data <- c(upstream = as.character(upstream),
                downstream = as.character(downstream),
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
get_flanking_seqs <- function(vr, up, down,
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
  mcols(vr) <- as.data.frame(do.call(rbind, flanking_seqs))

  return(vr)
}
