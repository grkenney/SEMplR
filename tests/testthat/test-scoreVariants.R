# set global variables for tests
b <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens

# define an SEM matrix to use in several tests
SEM_MATRIX <- t(matrix(c(0.000000, -0.7979680, -0.527404, -0.4609850, 
                               -0.288807, -0.6255560, -0.498335,  0.0000000,
                               -0.223680,  0.0795849, -0.401715, -0.0864365, 
                               0.000000, -0.5196070, -0.370991, -0.2774470,
                               0.000000, -0.4658680, -0.565092, -0.1388000, 
                               0.000000, -0.1555840, -0.482033,  0.1438920), 
                             nrow=4)) |> data.table()

colnames(SEM_MATRIX) <- c("A", "C", "G", "T")


test_that("scoreVariants SNP", {
  so <- SNPEffectMatrix(SEM_MATRIX, -0.402088,
                        "MA0151.1") 
  
  x <- VariantAnnotation::VRanges(seqnames = c("chr12", "chr19"),
                                  ranges = c(94136009, 10640062), 
                                  ref = c("G", "T"), alt = c("C", "A"))
  
  scores_a <- scoreVariants(x = x, sem = so, 
                            genome = b) |> 
    scores()
  
  scores_e <- data.table(varId = c('chr12:94136009:G>C', 'chr19:10640062:T>A'),
                         semId = c("MA0151.1"),
                         refSeq = c("TTTGAG", "ATCTCC"),
                         altSeq = c("TTCAGG", "ATCACC"),
                         refScore = c(-1.4004, -0.8193),
                         altScore = c(-1.4285, -0.5419),
                         refNorm = c(-0.4994, -0.2511),
                         altNorm = c(-0.5091, -0.0923),
                         refVarIndex = c(4, 4),
                         altVarIndex = c(5, 4))
  expect_equal(scores_a, scores_e, tolerance = 1e-3)
})


test_that("scoreVariants 1bp deletion", {
  so <- SNPEffectMatrix(SEM_MATRIX, -0.402088,
                        "MA0151.1")
  
  x <- VariantAnnotation::VRanges(seqnames = c("chr12"),
                                  ranges = c(94136009), 
                                  ref = c("G"), alt = c(""))
  
  scores_a <- scoreVariants(x = x, sem = so, 
                            genome = b) |> 
    scores()
  scores_a[, c("refScore", "altScore", "refNorm", "altNorm")] <-
    scores_a[, c("refScore", "altScore", "refNorm", "altNorm")] |>
    round(4)
  
  scores_e <- data.table(varId = c('chr12:94136009:delG'),
                         semId = c("MA0151.1"),
                         refSeq = c("TTTGAG"),
                         altSeq = c("TTTAGG"),
                         refScore = c(-1.4004),
                         altScore = c(-1.5945),
                         refNorm = c(-0.4994),
                         altNorm = c(-0.5624),
                         refVarIndex = c(4),
                         altVarIndex = c(4))
  expect_equal(scores_a, scores_e)
})


test_that("scoreVariants make semList a named list if not already", {
  so <- SNPEffectMatrix(SEM_MATRIX, -0.402088,
                        "MA0151.1")
  
  x <- VRanges(seqnames = c("chr12"),
                ranges = c(94136009), 
                ref = c("G"), alt = c(""))
  
  scores_a <- scoreVariants(x = x, sem = so, 
                            genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens) |> 
    scores()
  scores_a[, c("refScore", "altScore", "refNorm", "altNorm")] <-
    scores_a[, c("refScore", "altScore", "refNorm", "altNorm")] |>
    round(4)
  
  scores_e <- data.table(varId = c('chr12:94136009:delG'),
                         semId = c("MA0151.1"),
                         refSeq = c("TTTGAG"),
                         altSeq = c("TTTAGG"),
                         refScore = c(-1.4004),
                         altScore = c(-1.5945),
                         refNorm = c(-0.4994),
                         altNorm = c(-0.5624),
                         refVarIndex = c(4),
                         altVarIndex = c(4))
  expect_equal(scores_a, scores_e)
})


test_that("scoreVariants multiple variants not in alphanumeric order", {
  # create a matrix where the ideal sequence is ACGT
  SEM_MATRIX <- t(matrix(c(1, 0, 0, 0, 
                           0, 1, 0, 0,
                           0, 0, 1, 0,
                           0, 0, 0, 1), 
                         nrow=4)) |> data.table()
  colnames(SEM_MATRIX) <- c("A", "C", "G", "T")
  
  so <- SNPEffectMatrix(SEM_MATRIX, 1,
                        "sem_id")
  
  x <- VRanges(seqnames = c("chr12", "chr19"),
                ranges = c(94136009, 54282691), 
                ref = c("G", "G"), alt = c("", "T"))
  x$id <- c("B", "A")
  
  scores_a <- scoreVariants(x = x, sem = so, 
                            genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                            varId = "id") |> 
    scores()
  scores_a[, c("refScore", "altScore", "refNorm", "altNorm")] <-
    scores_a[, c("refScore", "altScore", "refNorm", "altNorm")] |>
    round(4)
  
  scores_e <- data.table(varId = c("A", "B"),
                         semId = c("sem_id"),
                         refSeq = c("ACGC", "TTGA"),
                         altSeq = c("ACTC", "TAGG"),
                         refScore = c(3, 1),
                         altScore = c(2, 1),
                         refNorm = c(3, 0),
                         altNorm = c(1, 0),
                         refVarIndex = c(3, 3),
                         altVarIndex = c(3, 4))
  expect_equal(scores_a, scores_e)
})


test_that("scoreVariants score a GRanges object", {
  so <- SNPEffectMatrix(SEM_MATRIX, -0.402088,
                        "MA0151.1")
  
  x <- GenomicRanges::GRanges(seqnames = c("chr12"),
                              ranges = c(94136009), 
                              ref = c("G"), alt = c(""))
  
  scores_a <- scoreVariants(x = x, sem = so, 
                            genome = b, 
                            refCol = "ref", altCol = "alt") |> 
    scores()
  scores_a[, c("refScore", "altScore", "refNorm", "altNorm")] <-
    scores_a[, c("refScore", "altScore", "refNorm", "altNorm")] |>
    round(4)
  
  scores_e <- data.table(varId = c('chr12:94136009:delG'),
                         semId = c("MA0151.1"),
                         refSeq = c("TTTGAG"),
                         altSeq = c("TTTAGG"),
                         refScore = c(-1.4004),
                         altScore = c(-1.5945),
                         refNorm = c(-0.4994),
                         altNorm = c(-0.5624),
                         refVarIndex = c(4),
                         altVarIndex = c(4))
  expect_equal(scores_a, scores_e)
})

