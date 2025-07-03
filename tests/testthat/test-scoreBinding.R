test_that(".makePositionId builds single position id", {
  x <- GenomicRanges::GRanges(seqnames = "chr12",
                              ranges = 94136009)
  id_a <- .makePositionId(x)
  expect_equal(id_a, "chr12:94136009")
})

test_that(".makePositionId builds single range id", {
  x <- GenomicRanges::GRanges(seqnames = "chr12",
                              ranges = IRanges::IRanges(94136009, 94136010))
  id_a <- .makePositionId(x)
  expect_equal(id_a, "chr12:94136009-94136010")
})

test_that(".makePositionId builds multiple ids", {
  x <- GenomicRanges::GRanges(seqnames = c("chr12", "chr13"),
                              ranges = IRanges::IRanges(c(94136009, 50000000), 
                                                        c(94136010, 50000000)))
  id_a <- .makePositionId(x)
  expect_equal(id_a, c("chr12:94136009-94136010", "chr13:50000000"))
})

# test_that(".scoreSingleSite single position for subset of SEMs", {
#   x <- GenomicRanges::GRanges(seqnames = "chr12",
#                               ranges = 94136009,
#                               seq = "CCGTCAAGGAGAAGGCTTTGAGGCATCTGCTGTTTTGTT")
#   s_a <- .scoreSingleSite(x = x, 
#                           semList = SNPEffectMatrixCollection(sems(sc)[1:2]), 
#                    bs_genome_obj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
#                    offset = 19)
#   s_e <- data.table::data.table(var_id = "chr12:94136009",
#                                 sem_mtx_id = c("AP2B_HUMAN.SK-N-SH",
#                                                "ARNT_HUMAN.GM12878"),
#                                 seq = c("GCTTTGAGGC", "TTTGAGGCA"),
#                                 vari = c(15, 17),
#                                 score = c(-1.689754, -6.892799),
#                                 scoreNorm = c(-0.30682393, -0.96938333))
#   expect_equal(s_a, s_e)
# })

test_that("scoreBinding 2 positions for subset of SEMs", {
  data(sc)
  x <- GenomicRanges::GRanges(seqnames = c("chr12", "chr13"),
                              ranges = IRanges::IRanges(c(94136009, 50000000), 
                                                        c(94136010, 50000000)),
                              seq = c("CCGTCAAGGAGAAGGCTTTGAGGCATCTGCTGTTTTGTT", 
                                      "GATCCATCTTCCAGGTTCAGCCCTCTGAGACACTCTTTC"))
  s_a <- scoreBinding(x = x, semList = sems(sc)[1:2], 
                      bs_genome_obj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
  
  # expcted data
  var_id_e <- rep(c("chr12:94136009-94136010", "chr13:50000000"), 2)
  sem_mtx_id_e <- rep(c("AP2B_HUMAN.SK-N-SH", "ARNT_HUMAN.GM12878"), each = 2)
  seq_e <- c("GCTTTGAGGC", "GGTTCAGCCC", "TTTGAGGCA", "GTTCAGCCC")
  score_e <- c(-1.689754, -2.983708, -6.892799, -6.997366)
  scoreNorm_e <- c(-0.30682393, -0.71730079, -0.96938333, -0.97152392)
  s_e <- data.table::data.table(var_id = var_id_e,
                                sem_mtx_id = sem_mtx_id_e,
                                seq = seq_e,
                                vari = c(6, 5, 8, 6),
                                score = score_e,
                                scoreNorm = scoreNorm_e)
  expect_equal(s_a, s_e)
})
