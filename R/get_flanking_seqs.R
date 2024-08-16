#' Get up and downstream sequences of a variant
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
#' @examples
#'
#' # Create a VRanges object
#' vr <- VRanges(seqnames = "chr1", ranges = 15000, ref = "A", alt = "T")
#'
#' # Get the 5 bp upstream and 10 bp downstream
#' get_flanking_seqs(vr, 5, 10)
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
  for (i in seq(vr)){
    # get start position of variant
    start_pos <- BiocGenerics::start(vr[i])
    end_pos <- BiocGenerics::end(vr[i])

    # get sequence upstream and downstream of variant start/end
    upstream <- Biostrings::getSeq(bs_genome_obj,
                       names=GenomeInfoDb::seqnames(vr[i]),
                       start=start_pos-up,
                       end=start_pos-1)
    downstream <- Biostrings::getSeq(bs_genome_obj,
                         names=GenomeInfoDb::seqnames(vr[i]),
                         start=end_pos+1,
                         end=end_pos+down)

    # store up and downstream sequences in vrange metadata
    vr[i]$upstream <- as.character(upstream)
    vr[i]$downstream <- as.character(downstream)

    # concatenate flanking sequences with variant and store in vrange metadata
    concat_ref_seq <- Biostrings::xscat(upstream, ref(vr[i]), downstream)
    concat_alt_seq <- Biostrings::xscat(upstream, alt(vr[i]), downstream)

    vr[i]$ref_seq <- as.character(concat_ref_seq)
    vr[i]$alt_seq <- as.character(concat_alt_seq)
  }
  return(vr)
}
