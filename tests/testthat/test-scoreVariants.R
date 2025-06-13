
# define an SEM matrix to use in several tests
SEM_MATRIX <- t(matrix(c(0.000000, -0.7979680, -0.527404, -0.4609850, 
                               -0.288807, -0.6255560, -0.498335,  0.0000000,
                               -0.223680,  0.0795849, -0.401715, -0.0864365, 
                               0.000000, -0.5196070, -0.370991, -0.2774470,
                               0.000000, -0.4658680, -0.565092, -0.1388000, 
                               0.000000, -0.1555840, -0.482033,  0.1438920), 
                             nrow=4)) |> data.table()

colnames(SEM_MATRIX) <- c("A", "C", "G", "T")


test_that("concatSeqs", {
  test_seq <- c("GCTTTGAGGCA",
                "GGATCTCCTGA")
  seq_concat_a <- concatSeqs(test_seq)
  seq_concat_e <- "GCTTTGAGGCAGGATCTCCTGA"
  expect_equal(seq_concat_a, seq_concat_e)
  
})


test_that("calcFrameStarts", {
  var_starts_a <- calcFrameStarts(nbp = 11, nflank = 5, motif_size = 3)
  expect_equal(var_starts_a, c(4, 5, 6))
})


test_that("normMatrix", {
  norm_a <- normMatrix(max_scores=c(13, 7), 
                       bl=-0.402088)
  norm_e <- c((2^13 - 2^-0.402088) / abs(2^-0.402088),
              (2^7 - 2^-0.402088) / abs(2^-0.402088))
  expect_equal(norm_e, norm_a)
  
})


test_that("scoreMatrix", {
  pwm <- data.table(A=c(1, 2), C=c(3, 4), G=c(5, 6), `T` = c(7, 8))
  
  scores_a <- scoreMatrix(sem_mtx = pwm, 
                          concat_seq = "ACTGTCAG",
                          frame_starts = list(c(2, 3), c(6, 7)))
  scores_e <- list(c(3+8, 7+6), c(3+2, 1+6))
  expect_equal(scores_a, scores_e)
})


test_that("calculateScores", {
  so <- SNPEffectMatrix(SEM_MATRIX, -0.402088,
                        "MA0151.1")
  scores_a <- calculateScores(varId = "rs13216361",
                              varSeq = "TGGCATTTCCTGGCAGAGCCCTGCCTCCCAGCTGTCTAA",
                              semObj = so,
                              nflank = 19)
  scores_e <- data.table(var_id = "rs13216361",
                         sem_mtx_id = "MA0151.1",
                         seq = c("AGCCCT"),
                         vari = c(17),
                         score = -1.2603331,
                         scoreNorm = -0.44837685)
  expect_equal(scores_a, scores_e)
})


test_that("scoreVariants SNP", {
  so <- SNPEffectMatrix(SEM_MATRIX, -0.402088,
                        "MA0151.1") 
  
  vr <- VRanges(seqnames = c("chr12", "chr19"),
                ranges = c(94136009, 10640062), 
                ref = c("G", "T"), alt = c("C", "A"))
  
  scores_a <- scoreVariants(vr, so, 
                            bs_genome_obj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens) |> 
    scores()
  scores_a[, c("refScore", "altScore", "refNorm", "altNorm")] <-
    scores_a[, c("refScore", "altScore", "refNorm", "altNorm")] |>
    round(4)
  
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
  expect_equal(scores_a, scores_e)
})


test_that("calcFrameStarts 1bp deletion", {
  fs_a <- calcFrameStarts(nbp = 4, nflank = 2, motif_size = 2)
  expect_equal(fs_a, 2)
})


test_that("calcFrameStarts 1bp insertion", {
  fs_a <- calcFrameStarts(nbp = 6, nflank = 2, motif_size = 2)
  expect_equal(fs_a, c(2, 3, 4))
})


test_that("calculateScores 1bp deletion", {
  so <- SNPEffectMatrix(SEM_MATRIX, -0.402088,
                        "MA0151.1")
  
  scores_a1 <- calculateScores("1", "GGCTTTGAGGCAT", so, nflank = 6)
  
  scores_e1 <- data.table(var_id = c('1'),
                          sem_mtx_id = c("MA0151.1"),
                          seq = c("TTTGAG"),
                          vari = c(4),
                          score = c(-1.4004455),
                          scoreNorm = c(-0.49943043))
  expect_equal(scores_a1, scores_e1)
  
  scores_a2 <- calculateScores("1", "GGCTTTAGGCAT", so, nflank = 6)
  
  scores_e2 <- data.table(var_id = c('1'),
                          sem_mtx_id = c("MA0151.1"),
                          seq = c("TTTAGG"),
                          vari = c(4),
                          score = c(-1.5945465),
                          scoreNorm = c(-0.56244342))
  expect_equal(scores_a2, scores_e2)
})


test_that("calculateScores multiple variants", {
  so <- SNPEffectMatrix(SEM_MATRIX, -0.402088,
                        "MA0151.1")
  
  scores_a1 <- calculateScores(c("1", "2", "3"), 
                             c("GGCTTTGAGGCAT", 
                               "GGCTTTAGGCAT", 
                               "GGCTTTGAGGCAT"), so, nflank = 6)
  
  scores_e1 <- data.table(var_id = c('1', '2', '3'),
                          sem_mtx_id = c("MA0151.1"),
                          seq = c("TTTGAG", "TTTAGG", "TTTGAG"),
                          vari = c(4),
                          score = c(-1.4004455, -1.5945465, -1.4004455),
                          scoreNorm = c(-0.49943043, -0.5624434, -0.49943043))
  expect_equal(scores_a1, scores_e1)
})


test_that("scoreVariants 1bp deletion", {
  so <- SNPEffectMatrix(SEM_MATRIX, -0.402088,
                        "MA0151.1")
  
  vr <- VRanges(seqnames = c("chr12"),
                ranges = c(94136009), 
                ref = c("G"), alt = c(""))
  
  scores_a <- scoreVariants(vr, so, 
                            bs_genome_obj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens) |> 
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


test_that("calculateScores variant length is longer than motif length", {
  so <- SNPEffectMatrix(SEM_MATRIX, -0.402088,
                        "MA0151.1")
  
  expect_error(calculateScores(varId = NA, varSeq = "ACTG", 
                               semObj = so, nflank = NA), 
               "Make sure all variant sequences are greater than or equal to")
})


test_that("scoreVariants make semList a named list if not already", {
  so <- SNPEffectMatrix(SEM_MATRIX, -0.402088,
                        "MA0151.1")
  
  vr <- VRanges(seqnames = c("chr12"),
                ranges = c(94136009), 
                ref = c("G"), alt = c(""))
  
  scores_a <- scoreVariants(vr, so, 
                            bs_genome_obj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens) |> 
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
  
  vr <- VRanges(seqnames = c("chr12", "chr19"),
                ranges = c(94136009, 54282691), 
                ref = c("G", "G"), alt = c("", "C"))
  vr$id <- c("B", "A")
  
  scores_a <- scoreVariants(vr, so, 
                            bs_genome_obj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens) |> 
    scores()
  scores_a[, c("refScore", "altScore", "refNorm", "altNorm")] <-
    scores_a[, c("refScore", "altScore", "refNorm", "altNorm")] |>
    round(4)
  
  scores_e <- data.table(varId = c("A", "B"),
                         semId = c("sem_id"),
                         refSeq = c("ACGC", "TTGA"),
                         altSeq = c("ACCC", "TAGG"),
                         refScore = c(3, 1),
                         altScore = c(2, 1),
                         refNorm = c(3, 0),
                         altNorm = c(1, 0),
                         refVarIndex = c(3, 3),
                         altVarIndex = c(3, 4 ))
  expect_equal(scores_a, scores_e)
})
