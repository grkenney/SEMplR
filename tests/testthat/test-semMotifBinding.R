
# define an SEM matrix to use in several tests
SEM_MATRIX <- t(matrix(c(0.000000, -0.7979680, -0.527404, -0.4609850, 
                               -0.288807, -0.6255560, -0.498335,  0.0000000,
                               -0.223680,  0.0795849, -0.401715, -0.0864365, 
                               0.000000, -0.5196070, -0.370991, -0.2774470,
                               0.000000, -0.4658680, -0.565092, -0.1388000, 
                               0.000000, -0.1555840, -0.482033,  0.1438920), 
                             nrow=4)) |> data.table()

colnames(SEM_MATRIX) <- c("A", "C", "G", "T")

test_that("concatVars", {
  test_seq <- c("GCTTT-G-AGGCA",
                "GGATC-T-CCTGA")
  seq_concat_a <- concatVars(test_seq)
  seq_concat_e <- "GCTTT-G-AGGCAGGATC-T-CCTGA"
  expect_equal(seq_concat_a, seq_concat_e)
  
})

test_that("calcVarStarts", {
  var_starts_a <- calcVarStarts(offset = 5, 
                                nVars = 2, 
                                nFrames = 3)
  expect_equal(var_starts_a, c(4, 5, 6, 15, 16, 17))
})

test_that("normMatrix", {
  norm_a <- normMatrix(scores_mtx=matrix(c(11, 13, 5, 7), 
                                         nrow=2), 
                       bl=-0.402088)
  norm_e <- c((2^13 - 2^-0.402088) / abs(2^-0.402088),
              (2^7 - 2^-0.402088) / abs(2^-0.402088))
  expect_equal(norm_e, norm_a)
  
})

test_that("scoreMatrix", {
  pwm <- data.table(A=c(1, 2), C=c(3, 4), G=c(5, 6), `T` = c(7, 8))
  
  scores_a <- scoreMatrix(sem_mtx = pwm, 
                          concat_seq = "ACTGTCAG",
                          frame_starts = c(2, 3, 6, 7),
                          nFrames = 2)
  scores_e <- matrix(c(3+8, 7+6, 3+2, 1+6), nrow=2)
  expect_equal(scores_a, scores_e)
})


test_that("scoreVariants", {
  so <- SNPEffectMatrix(SEM_MATRIX, -0.402088,
                        "MA0151.1", tf = "ARID3A", 
                        ensembl = "ENSG00000116017", uniprot = "Q99856",
                        cellType = "HepG2")
  scores_a <- scoreVariants(varId = "rs13216361",
                            varSeq = "TGGCATTTCCTGGCAGAGCCCTGCCTCCCAGCTGTCTAA",
                            semObj = so,
                            offset = 19)
  scores_e <- data.table(var_id = "rs13216361",
                         sem_mtx_id = "MA0151.1",
                         seq = "AGCCCT",
                         vari = 4,
                         score = -1.2603331,
                         scoreNorm = -0.44837685)
  expect_equal(scores_a, scores_e)
})
