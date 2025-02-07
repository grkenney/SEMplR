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