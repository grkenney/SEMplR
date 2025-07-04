test_that("querySeqs pulls correct seqs", {
  vr <- VRanges(seqnames = "chr12",
                ranges = 94136009,
                ref = "G", alt = "C")
  seq_data_e <- c(upstream = "TT",
                  downstream = "AG",
                  ref_seq = "TTGAG",
                  alt_seq = "TTCAG")
  seq_data_a <- querySeqs(vr, up = 2, down = 2,
                          bs_genome_obj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
  expect_equal(seq_data_a, seq_data_e)
})


test_that("getFlankingSeqs normal vranges input", {
  vr <- VRanges(seqnames = "chr12",
                ranges = 94136009,
                ref = "G", alt = "C")
  vr_a <- getFlankingSeqs(vr, 5, 5, 
                          bs_genome_obj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
  
  # populate VRanges and with expected sequences
  vr$upstream <- "GCTTT"
  vr$downstream <- "AGGCA"
  vr$ref_seq <- "GCTTTGAGGCA"
  vr$alt_seq <- "GCTTTCAGGCA"
  
  expect_equal(vr_a, vr)
})


test_that("getFlankingSeqs vr not of class VRanges", {
  expect_error(getFlankingSeqs(NA, 5, 5, 
                               bs_genome_obj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens), 
               "Input argument x must be of class VRanges or GRanges")
})


test_that("getFlankingSeqs VRanges is empty", {
  expect_error(getFlankingSeqs(VRanges(), 5, 5), 
               "VRanges object must contain at least one")
})
