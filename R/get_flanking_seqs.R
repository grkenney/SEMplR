library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

# get up and downstream seqs
get_flanking_seqs <- function(vr, up, down, bs_genome_obj=Hsapiens){
  # initialize metadata columns in vrange to store up and downstream seqs
  vr$upstream <- NA
  vr$downstream <- NA

  # iterate over variants (rows) in vrange object
  for (i in seq(vr)){
    # get start position of variant
    start_pos <- start(vr[i])
    end_pos <- end(vr[i])

    # get sequence upstream and downstream of variant start/end
    upstream <- getSeq(bs_genome_obj,
                       names=seqnames(vr[i]),
                       start=start_pos-up,
                       end=start_pos-1)
    downstream <- getSeq(bs_genome_obj,
                         names=seqnames(vr[i]),
                         start=end_pos+1,
                         end=end_pos+down)

    # store up and downstream sequences in vrange metadata
    vr[i]$upstream <- as.character(upstream)
    vr[i]$downstream <- as.character(downstream)
  }
  return(vr)
}
