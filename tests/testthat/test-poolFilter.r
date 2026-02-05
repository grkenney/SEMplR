orgdb <- org.Hs.eg.db::org.Hs.eg.db
my_genes <- c("ENSG00000139618", "ENSG00000157764")


test_that("poolFilter minimal example", {
  ids <- mapIDs(
    orgdb = orgdb,
    foreground_ids = my_genes,
    id_type = "ENSEMBL")
  pf_a <- poolFilter(ids, geneType = "protein-coding")
  pf_e <- data.frame(ENSEMBL = my_genes, ENTREZID = c("675", "673"), 
                         GENETYPE = "protein-coding")
  expect_equal(pf_a$fg_ids, pf_e)
})


test_that("poolFilter invalid input", {
  ids <- mapIDs(
    orgdb = orgdb,
    foreground_ids = my_genes,
    id_type = "ENSEMBL")
  expect_error(poolFilter(ids, geneType = "foo"), 
               regexp = "Invalid geneType 'foo'")
})
