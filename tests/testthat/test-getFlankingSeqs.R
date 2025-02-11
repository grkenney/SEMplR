test_that("query_seqs pulls correct seqs", {
  vr <- VRanges(seqnames = "chr12",
                ranges = 94136009,
                ref = "G", alt = "C")
  seq_data_e <- c(upstream = "TT",
                  downstream = "AG",
                  ref_seq = "TTGAG",
                  alt_seq = "TTCAG")
  seq_data_a <- query_seqs(vr, up = 2, down = 2)
  expect_equal(seq_data_a, seq_data_e)
})


test_that("getFlankingSeqs normal vranges input", {
  vr <- VRanges(seqnames = "chr12",
                ranges = 94136009,
                ref = "G", alt = "C")
  vr_a <- getFlankingSeqs(vr, 5, 5)
  
  # populate VRanges and with expected sequences
  vr$upstream <- "GCTTT"
  vr$downstream <- "AGGCA"
  vr$ref_seq <- "GCTTTGAGGCA"
  vr$alt_seq <- "GCTTTCAGGCA"
  
  expect_equal(vr_a, vr)
})


test_that("getFlankingSeqs vr not of class VRanges", {
  expect_error(getFlankingSeqs(NA, 5, 5), 
               "Input argument vr must be of class VRanges")
})


test_that("getFlankingSeqs VRanges is empty", {
  expect_error(getFlankingSeqs(VRanges(), 5, 5), 
               "VRanges object must contain at least one")
})
