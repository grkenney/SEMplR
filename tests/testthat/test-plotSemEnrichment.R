test_that("prepPlotSemEnrichmentDf significant", {
  e <- data.table(semId = "sem_id",
                  tf = "tf_name",
                  n.changed = 1,
                  odds.ratio = 2,
                  ci.lower = 1,
                  ci.upper = 1,
                  pvalue = 0.05,
                  adj.pvalue = 0.05)
  e_a <- prepPlotSemEnrichmentDf(e, lab = "tf", sigThreshold = 0.05)
  e$lab <- "tf_name"
  e$sig <- factor("adj. p-value ≤ 0.05", 
                  levels = c("adj. p-value ≤ 0.05", "> 0.05"))
  expect_equal(e_a, e)
})


test_that("prepPlotSemEnrichmentDf not significant", {
  e <- data.table(semId = "sem_id",
                  tf = "tf_name",
                  n.changed = 1,
                  odds.ratio = 1,
                  ci.lower = 1,
                  ci.upper = 1,
                  pvalue = 0.05,
                  adj.pvalue = 0.1)
  e_a <- prepPlotSemEnrichmentDf(e, lab = "tf", sigThreshold = 0.05)
  e$lab <- ""
  e$sig <- factor("> 0.05", 
                  levels = c("adj. p-value ≤ 0.05", "> 0.05"))
  expect_equal(e_a, e)
})
