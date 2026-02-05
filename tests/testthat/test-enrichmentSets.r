orgdb <- org.Hs.eg.db::org.Hs.eg.db
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
my_genes <- c("ENSG00000139618", "ENSG00000157764")
bg_genes <- c("ENSG00000109321", "ENSG00000078061")


test_that("enrichmentSets minimal example", {
  es <- enrichmentSets(txdb = txdb,
                       orgdb = orgdb,
                       id_type = "ENSEMBL",
                       foreground_ids = my_genes,
                       background_ids = bg_genes)
  expect_in(as.vector(GenomicRanges::seqnames(es$foregroundElements)), 
            c("chr13", "chr7"))
  expect_in(GenomicRanges::start(es$foregroundElements), 
            c(32314786, 140924880))
  expect_in(as.vector(GenomicRanges::seqnames(es$backgroundElements)), 
            c("chrX", "chr4"))
  expect_in(GenomicRanges::start(es$backgroundElements), 
            c(47560905, 74444836))
})


test_that("mapIds errors and messages", {
  expect_error(mapIDs(
    orgdb = orgdb,
    foreground_ids = my_genes,
    id_type = "foo"), 
    regexp = "not an available key for in this mapping object")
  expect_message(mapIDs(
    orgdb = orgdb,
    foreground_ids = my_genes,
    background_ids = "ENSG00000139610",
    id_type = "ENSEMBL"), 
    regexp = "Mapping background ids")
})
