#' Builds a VRanges object from genomic position and allele information
#'
#' @param chr character, chromosome identifier
#' @param pos integer, SNP position on chromosome
#' @param ref character, reference nucleotide (can be upper or lower)
#' @param alt character, alternate nulceotide (can be upper or lower)
#'
#' @return a `VRanges` object the seqnames, ranges, ref, and alt columns
#' filled with parameter options supplied.
#'
#' @examples
#'
#' vr <- VRanges(seqnames = "chr1", ranges = 15000, ref = "A", alt = "T")
#' vr
#'
#' @export
chr_pos_to_vrange <- function(chr, pos, ref, alt){
  vr <- VRanges(chr, ranges = pos,
                ref = ref, alt = alt)
  return(vr)
}
