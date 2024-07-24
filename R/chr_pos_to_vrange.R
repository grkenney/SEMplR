library(GenomicRanges)

# function that makes vranges from genomic position
chr_pos_to_vrange <- function(chr, pos, ref, alt){
  vr <- VRanges(chr, ranges = pos,
                ref = ref, alt = alt)
  return(vr)
}
