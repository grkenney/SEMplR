SEM_MATRIX <- t(matrix(1:24, nrow=4)) |> 
  data.table()
colnames(SEM_MATRIX) <- c("A", "C", "G", "T")

SEM_OBJ <- SNPEffectMatrix(SEM_MATRIX, -0.402088,
                           "MA0151.1") |> 
  list() |>
  setNames("MA0151.1")

VR <- VRanges(seqnames = "1",
              ranges = IRanges::IRanges(start = 1:8, end = 1:8),
              ref = c("A", "", "A", "AA", "", "AA", "A", "A"), 
              alt = c("T", "T", "", "", "TT", "TT", "T", "T"))
VR$id <- 1:8

VR$upstream <- "CCCCCC"
VR$downstream <- "GGGGGG"

SCORES_DT <- data.table(varId=1:8, 
                        semId="MA0151.1", 
                        refSeq="CCCAGG", 
                        altSeq="CCCTGG",
                        refScore=NA, 
                        altScore=NA,
                        refNorm=NA, 
                        altNorm=NA,
                        refVarIndex=c(3:9, 8),
                        altVarIndex=c(3:9, 8))

SEMPL_SCORES_OBJ <- SEMplScores(ranges = VR,
                                semData = data.table(),
                                scores = SCORES_DT)


test_that("viewFrames 1bp SNP", {
  vf_a <- viewFrames(s = SEMPL_SCORES_OBJ, scoreIndex = 1)
  vf_e <- SequenceFrame(seqName = c("ref", "alt"),
                        frameStart = c(3, 3), 
                        frameStop = c(8, 8), 
                        sequence = c("CCCCCCAGGGGGG", "CCCCCCTGGGGGG"),
                        variantIndex = 7,
                        motif = "MA0151.1",
                        variantName = "1")
  expect_equal(vf_a, vf_e)
})


test_that("viewFrames 1bp insertion", {
  vf_a <- viewFrames(s = SEMPL_SCORES_OBJ, scoreIndex = 2)
  vf_e <- SequenceFrame(seqName = c("ref", "alt"),
                        frameStart = c(4, 4), 
                        frameStop = c(9, 9), 
                        sequence = c("CCCCCC-GGGGGG", "CCCCCCTGGGGGG"),
                        variantIndex = 7,
                        motif = "MA0151.1",
                        variantName = "2")
  expect_equal(vf_a, vf_e)
})


test_that("viewFrames 1bp deletion", {
  vf_a <- viewFrames(s = SEMPL_SCORES_OBJ, scoreIndex = 3)
  vf_e <- SequenceFrame(seqName = c("ref", "alt"),
                        frameStart = c(5, 5), 
                        frameStop = c(10, 10), 
                        sequence = c("CCCCCCAGGGGGG", "CCCCCC-GGGGGG"),
                        variantIndex = 7,
                        motif = "MA0151.1",
                        variantName = "3")
  expect_equal(vf_a, vf_e)
})


test_that("viewFrames 2bp deletion", {
  vf_a <- viewFrames(s = SEMPL_SCORES_OBJ, scoreIndex = 4)
  vf_e <- SequenceFrame(seqName = c("ref", "alt"),
                        frameStart = c(6, 6), 
                        frameStop = c(11, 11), 
                        sequence = c("CCCCCCAAGGGGGG", "CCCCCC--GGGGGG"),
                        variantIndex = 7:8,
                        motif = "MA0151.1",
                        variantName = "4")
  expect_equal(vf_a, vf_e)
})


test_that("viewFrames 2bp insertion", {
  vf_a <- viewFrames(s = SEMPL_SCORES_OBJ, scoreIndex = 5)
  vf_e <- SequenceFrame(seqName = c("ref", "alt"),
                        frameStart = c(7, 7), 
                        frameStop = c(12, 12), 
                        sequence = c("CCCCCC--GGGGGG", "CCCCCCTTGGGGGG"),
                        variantIndex = 7:8,
                        motif = "MA0151.1",
                        variantName = "5")
  expect_equal(vf_a, vf_e)
})


test_that("viewFrames 2bp SNP", {
  vf_a <- viewFrames(s = SEMPL_SCORES_OBJ, scoreIndex = 6)
  vf_e <- SequenceFrame(seqName = c("ref", "alt"),
                        frameStart = c(8, 8), 
                        frameStop = c(13, 13), 
                        sequence = c("CCCCCCAAGGGGGG", "CCCCCCTTGGGGGG"),
                        variantIndex = 7:8,
                        motif = "MA0151.1",
                        variantName = "6")
  expect_equal(vf_a, vf_e)
})


test_that("viewFrames SNP not in frame", {
  expect_error(viewFrames(s = SEMPL_SCORES_OBJ, scoreIndex = 7),
               "Frame index range")
})


test_that("viewFrames variant not in frame", {
  expect_warning(viewFrames(s = SEMPL_SCORES_OBJ, scoreIndex = 8),
                 "Variant index")
})


test_that("viewFrames scoreIndex not in scores", {
  expect_error(viewFrames(s = SEMPL_SCORES_OBJ, scoreIndex = 9),
               "out of range in SEMplScores object")
})


test_that("viewFrames no entry matching vid and sid", {
  expect_error(viewFrames(s = SEMPL_SCORES_OBJ, vid = 1, sid = 1),
               "no entry matching")
})


test_that("viewFrames require sid and vid or scoreIndex - neither", {
  # neither
  expect_error(viewFrames(s = SEMPL_SCORES_OBJ), 
               "must provide either a scoreIndex or a vid and sid")
})


test_that("viewFrames require sid and vid or scoreIndex - vid only", {
  # just vid
  expect_error(viewFrames(s = SEMPL_SCORES_OBJ, vid = 1), 
               "must provide either a scoreIndex or a vid and sid")
})


test_that("viewFrames require sid and vid or scoreIndex - sid only", {
  # just sid
  expect_error(viewFrames(s = SEMPL_SCORES_OBJ, sid = 1), 
               "must provide either a scoreIndex or a vid and sid")
})
