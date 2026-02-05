txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
orgdb <- org.Hs.eg.db::org.Hs.eg.db
my_genes <- c("ENSG00000139618", "ENSG00000157764")


test_that("getCoordinates minimal example", {
  ids <- mapIDs(
    orgdb = orgdb,
    foreground_ids = my_genes,
    background_ids = c("ENSG00000109321", "ENSG00000078061"),
    id_type = "ENSEMBL")
  pf <- poolFilter(ids, geneType = "protein-coding")
  coords_a <- getCoordinates(mapped = pf, txdb = txdb, n_ratio = 1)
  expect_in(as.vector(GenomicRanges::seqnames(coords_a$foregroundElements)), 
               c("chr13", "chr7"))
  expect_in(GenomicRanges::start(coords_a$foregroundElements), 
               c(32314786, 140924880))
  expect_in(as.vector(GenomicRanges::seqnames(coords_a$backgroundElements)), 
               c("chrX", "chr4"))
  expect_in(GenomicRanges::start(coords_a$backgroundElements), 
               c(47560905, 74444836))
})


test_that("getCoordinates invalid input", {
  ids <- mapIDs(
    orgdb = orgdb,
    foreground_ids = my_genes,
    background_ids = "ENSG00000139610",
    id_type = "ENSEMBL")
  pf <- poolFilter(ids, geneType = "protein-coding")
  expect_error(getCoordinates(mapped = pf, txdb = txdb),
               regexp = "Requested 2 background regions")
})

