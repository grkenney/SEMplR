orgdb <- org.Hs.eg.db::org.Hs.eg.db
my_genes <- c("ENSG00000139618", "ENSG00000157764")


test_that("mapIds minimal example", {
  ids_a <- mapIDs(
    orgdb = orgdb,
    foreground_ids = my_genes,
    id_type = "ENSEMBL")
  fg_ids_e <- data.frame(ENSEMBL = my_genes, ENTREZID = c("675", "673"), 
                         GENETYPE = "protein-coding")
  expect_equal(ids_a$fg_ids, fg_ids_e)
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
