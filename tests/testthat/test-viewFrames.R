SEM_MATRIX <- t(matrix(1:24, nrow=4)) |> 
  data.table()
colnames(SEM_MATRIX) <- c("A", "C", "G", "T")

SEM_OBJ <- SNPEffectMatrix(SEM_MATRIX, -0.402088,
                           "MA0151.1", tf = "ARID3A", 
                           ensembl = "ENSG00000116017", uniprot = "Q99856",
                           cellType = "HepG2") |> 
  list() |>
  setNames("MA0151.1")

VR <- VRanges(seqnames = "1",
              ranges = IRanges::IRanges(1),
              ref = c("A", "", "A", "AA", "", "AA", "A", "A"), 
              alt = c("T", "T", "", "", "TT", "TT", "T", "T"))

VR$upstream <- "CCCCCC"
VR$downstream <- "GGGGGG"

SCORES_DT <- data.table(varId=1:8, 
                        semId="MA0151.1", 
                        nonRiskSeq="CCCAGG", 
                        riskSeq="CCCTGG",
                        nonRiskScore=NA, 
                        riskScore=NA,
                        nonRiskNorm=NA, 
                        riskNorm=NA,
                        nonRiskVarIndex=c(3:9, 8),
                        riskVarIndex=c(3:9, 8))

SEMPL_SCORES_OBJ <- SemplScores(variants = VR,
                                semList = SEM_OBJ,
                                semScores = SCORES_DT)


test_that("viewFrames 1bp SNP", {
  vf_a <- viewFrames(sempl_obj = SEMPL_SCORES_OBJ, score_index = 1)
  vf_e <- SequenceFrame(seqName = c("ref", "alt"),
                        frameStart = c(3, 3), 
                        frameStop = c(8, 8), 
                        sequence = c("CCCCCCAGGGGGG", "CCCCCCTGGGGGG"),
                        variantIndex = 7,
                        motif = "MA0151.1",
                        tf = "ARID3A",
                        variantName = "1")
  expect_equal(vf_a, vf_e)
})


test_that("viewFrames 1bp insertion", {
  vf_a <- viewFrames(sempl_obj = SEMPL_SCORES_OBJ, score_index = 2)
  vf_e <- SequenceFrame(seqName = c("ref", "alt"),
                        frameStart = c(4, 4), 
                        frameStop = c(9, 9), 
                        sequence = c("CCCCCC-GGGGGG", "CCCCCCTGGGGGG"),
                        variantIndex = 7,
                        motif = "MA0151.1",
                        tf = "ARID3A",
                        variantName = "2")
  expect_equal(vf_a, vf_e)
})


test_that("viewFrames 1bp deletion", {
  vf_a <- viewFrames(sempl_obj = SEMPL_SCORES_OBJ, score_index = 3)
  vf_e <- SequenceFrame(seqName = c("ref", "alt"),
                        frameStart = c(5, 5), 
                        frameStop = c(10, 10), 
                        sequence = c("CCCCCCAGGGGGG", "CCCCCC-GGGGGG"),
                        variantIndex = 7,
                        motif = "MA0151.1",
                        tf = "ARID3A",
                        variantName = "3")
  expect_equal(vf_a, vf_e)
})


test_that("viewFrames 2bp deletion", {
  vf_a <- viewFrames(sempl_obj = SEMPL_SCORES_OBJ, score_index = 4)
  vf_e <- SequenceFrame(seqName = c("ref", "alt"),
                        frameStart = c(6, 6), 
                        frameStop = c(11, 11), 
                        sequence = c("CCCCCCAAGGGGGG", "CCCCCC--GGGGGG"),
                        variantIndex = 7:8,
                        motif = "MA0151.1",
                        tf = "ARID3A",
                        variantName = "4")
  expect_equal(vf_a, vf_e)
})


test_that("viewFrames 2bp insertion", {
  vf_a <- viewFrames(sempl_obj = SEMPL_SCORES_OBJ, score_index = 5)
  vf_e <- SequenceFrame(seqName = c("ref", "alt"),
                        frameStart = c(7, 7), 
                        frameStop = c(12, 12), 
                        sequence = c("CCCCCC--GGGGGG", "CCCCCCTTGGGGGG"),
                        variantIndex = 7:8,
                        motif = "MA0151.1",
                        tf = "ARID3A",
                        variantName = "5")
  expect_equal(vf_a, vf_e)
})


test_that("viewFrames 2bp SNP", {
  vf_a <- viewFrames(sempl_obj = SEMPL_SCORES_OBJ, score_index = 6)
  vf_e <- SequenceFrame(seqName = c("ref", "alt"),
                        frameStart = c(8, 8), 
                        frameStop = c(13, 13), 
                        sequence = c("CCCCCCAAGGGGGG", "CCCCCCTTGGGGGG"),
                        variantIndex = 7:8,
                        motif = "MA0151.1",
                        tf = "ARID3A",
                        variantName = "6")
  expect_equal(vf_a, vf_e)
})


test_that("viewFrames SNP not in frame", {
  expect_error(viewFrames(sempl_obj = SEMPL_SCORES_OBJ, score_index = 7),
               "Frame index range")
})


test_that("viewFrames variant not in frame", {
  expect_warning(viewFrames(sempl_obj = SEMPL_SCORES_OBJ, score_index = 8),
                 "Variant index")
})


test_that("viewFrames score_index not in scores", {
  expect_error(viewFrames(sempl_obj = SEMPL_SCORES_OBJ, score_index = 9),
               "out of range in SemplScores object")
})


test_that("viewFrames no entry matching vid and sid", {
  expect_error(viewFrames(sempl_obj = SEMPL_SCORES_OBJ, vid = 1, sid = 1),
               "no entry matching")
})


test_that("viewFrames require sid and vid or score_index - neither", {
  # neither
  expect_error(viewFrames(sempl_obj = SEMPL_SCORES_OBJ), 
               "must provide either a score_index or a vid and sid")
})


test_that("viewFrames require sid and vid or score_index - vid only", {
  # just vid
  expect_error(viewFrames(sempl_obj = SEMPL_SCORES_OBJ, vid = 1), 
               "must provide either a score_index or a vid and sid")
})


test_that("viewFrames require sid and vid or score_index - sid only", {
  # just sid
  expect_error(viewFrames(sempl_obj = SEMPL_SCORES_OBJ, sid = 1), 
               "must provide either a score_index or a vid and sid")
})
