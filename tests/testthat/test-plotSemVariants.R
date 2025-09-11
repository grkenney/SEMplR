# setup test example
x <- VRanges(seqnames = "chr12",
             ranges = 94136009, 
             ref = "G", alt = "C",
             id = "1")
s <- scoreVariants(x = x, sem = sc, 
                   genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens, 
                   varId = "id")

test_that(".validatePlotSemVariantsInputs errors on invalid input", {
  # no error
  expect_no_error(.validatePlotSemVariantsInputs(s = s, 
                                                 label = "id", 
                                                 semId = "MA0099.2_HeLa",
                                                 cols = c("red", "blue")))
  
  # invalid label
  expect_error(.validatePlotSemVariantsInputs(s = s, 
                                              label = "foo", 
                                              semId = "MA0099.2_HeLa",
                                              cols = c("red", "blue")),
               regexp = "label not found")
  
  # invalid semId
  expect_error(.validatePlotSemVariantsInputs(s = s, 
                                              label = "id", 
                                              semId = "foo",
                                              cols = c("red", "blue")),
               regexp = "variant not found")
  
  # invalid cols
  expect_error(.validatePlotSemVariantsInputs(s = s, 
                                              label = "id", 
                                              semId = "MA0099.2_HeLa",
                                              cols = c("red")),
               regexp = "cols must be a vector of length 2")
  expect_error(.validatePlotSemVariantsInputs(s = s, 
                                              label = "id", 
                                              semId = "MA0099.2_HeLa",
                                              cols = c("red", "blue", "green")),
               regexp = "cols must be a vector of length 2")
  
})


test_that("plotSemVariants returns a ggplot object", {
  plt <- plotSemVariants(s = s, 
                         sem = "AP2B_HUMAN.SK-N-SH")
  expect_s3_class(plt, "gg")
})
