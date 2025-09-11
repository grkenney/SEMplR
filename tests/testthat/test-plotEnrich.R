# setup enrichment
x <- GenomicRanges::GRanges(seqnames = c("chr12", "chr16"),
                            ranges = IRanges::IRanges(start = c(94136009,
                                                                67637901),
                                                      end = c(94136009,
                                                              67637901)))
b <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
sb <- scoreBinding(x = x, sem = sc, genome = b)
e <- enrichSEMs(x = sb, sem = sc, genome = b)


test_that(".formatMotifs adds correct altname", {
  motifs <- .formatMotifs(sem = sc, 
                          label = "transcription_factor")
  expect_equal(motifs[[1]]["altname"], "TFAP2B")
  expect_equal(motifs[[1]]["name"], "AP2B_HUMAN.SK-N-SH")
  expect_equal(motifs[[100]]["altname"], "GATA1")
  expect_equal(motifs[[100]]["name"], "M00347")
})


test_that(".constructComparisons returns correctly sized object", {
  motifs <- .formatMotifs(sem = sc, 
                          label = "transcription_factor")
  labels <- lapply(motifs, function(x) x["altname"]) |> unlist()
  
  comps <- .constructComparisons(motifs = motifs, labels = labels, 
                                 method = "WPCC")
  expect_s3_class(comps, "phylo")
  expect_length(comps$edge, 888)
})


test_that(".plotMotifTree returns tree object", {
  motifs <- .formatMotifs(sem = sc, 
                          label = "transcription_factor")
  labels <- lapply(motifs, function(x) x["altname"]) |> unlist()
  
  comps <- .constructComparisons(motifs = motifs, labels = labels, 
                                 method = "PCC")
  
  plt <- .plotMotifTree(e = e, sem = sc, threshold = 0.05, 
                        comparisons = comps, sigCol = "purple4")
  expect_s3_class(plt, "ggtree")
})


test_that(".plotTreeHeatmap returns ggtree object", {
  motifs <- .formatMotifs(sem = sc, 
                          label = "transcription_factor")
  labels <- lapply(motifs, function(x) x["altname"]) |> unlist()
  
  comps <- .constructComparisons(motifs = motifs, 
                                 labels = labels, 
                                 method = "PCC")
  
  plt <- .plotMotifTree(e = e, sem = sc, threshold = 0.05, 
                        comparisons = comps, sigCol = "purple4")
  plt <- .plotTreeHeatmap(plt = plt, e = e, labels = labels)
  expect_s3_class(plt, "ggtree")
})


test_that("plotEnrich generates a ggplot object", {
  plt <- plotEnrich(e, sc)
  expect_s3_class(plt, "ggtree")
})

