test_that(".validateMotifInCollection validates motifs correctly", {
    expect_no_error(.validateMotifInCollection(SEMC, "MA0099.2_HeLa"))
    expect_error(.validateMotifInCollection(SEMC, NULL),
        regexp = "must specify the motif"
    )
    expect_error(.validateMotifInCollection(SEMC, "foo"),
        regexp = "provided motif not found"
    )
})


test_that(".definePlotSEMParams works for SNPEffectMatrices", {
    sem_obj <- getSEMs(SEMC, "MA0099.2_HeLa")
    params_sem_a <- .definePlotSEMParams(
        sem = sem_obj,
        motif = NULL
    ) # without motif
    expect_equal(params_sem_a$sem_baseline, getBaseline(sem_obj))

    params_sem_a <- .definePlotSEMParams(
        sem = sem_obj,
        motif = "MA0099.2"
    ) # with motif
    expect_equal(params_sem_a$sem_baseline, getBaseline(sem_obj))
    expect_equal(params_sem_a$sem_mtx, getSEM(sem_obj))
    expect_message(.definePlotSEMParams(sem = sem_obj, motif = "MA0099.2"),
        regexp = "motif ignored"
    )
})


test_that(".definePlotSEMParams works for SNPEffectMatrixCollections", {
    sem_obj <- getSEMs(SEMC, "MA0099.2_HeLa")
    params_sem_a <- .definePlotSEMParams(sem = SEMC, motif = "MA0099.2_HeLa")
    expect_equal(params_sem_a$sem_baseline, getBaseline(sem_obj))
    expect_equal(params_sem_a$sem_mtx, getSEM(sem_obj))
})


test_that(".definePlotSEMParams fails on unexpected sem class", {
    # sem is not a SNPEffectMatrix of a Collection
    expect_error(.definePlotSEMParams(sem = TRUE, motif = NULL),
        regexp = "sem must be a SNPEffectMatrix or a"
    )
    expect_error(.definePlotSEMParams(sem = NULL, motif = NULL),
        regexp = "sem must be a SNPEffectMatrix or a"
    )
    expect_error(.definePlotSEMParams(sem = "foo", motif = NULL),
        regexp = "sem must be a SNPEffectMatrix or a"
    )
})


test_that(".formatPlotSEMTable pivots SEM data", {
    sem_mtx <- getSEMs(SEMC, "MA0099.2_HeLa") |>
        getSEM()
    sem_mtx_long <- .formatPlotSEMTable(sem_mtx)
    # has the correct number of rows
    expect_equal(nrow(sem_mtx_long), nrow(sem_mtx) * ncol(sem_mtx))
    # column names are correct
    expect_equal(colnames(sem_mtx_long), c("sem_score", "motif_pos", "bp"))
})


test_that(".defineNucleotideColors chooses black text when motifSeq is NULL", {
    sem_mtx <- getSEMs(SEMC, "MA0099.2_HeLa") |>
        getSEM()
    sem_mtx_long <- .formatPlotSEMTable(sem_mtx)

    res_a <- .defineNucleotideColors(
        sem_mtx_long = sem_mtx_long,
        motifSeq = NULL
    )
    expect_equal(colnames(res_a), c(
        "sem_score", "motif_pos",
        "bp", "text_color"
    ))
    expect_equal(res_a$text_color, rep("black", nrow(sem_mtx_long)))
})


test_that(".defineNucleotideColors fails when incorrect number of bases", {
    sem_mtx <- getSEMs(SEMC, "MA0099.2_HeLa") |>
        getSEM()
    sem_mtx_long <- .formatPlotSEMTable(sem_mtx)

    # too long
    expect_error(
        .defineNucleotideColors(
            sem_mtx_long = sem_mtx_long,
            motifSeq = "AAAAAAAA"
        ),
        regexp = "Number of characters in motifSeq"
    )
    # too short
    expect_error(
        .defineNucleotideColors(
            sem_mtx_long = sem_mtx_long,
            motifSeq = "AAAAAA"
        ),
        regexp = "Number of characters in motifSeq"
    )
})


test_that(".defineNucleotideColors gives text_colors for motifSeq", {
    sem_mtx <- getSEMs(SEMC, "MA0099.2_HeLa") |>
        getSEM()
    sem_mtx_long <- .formatPlotSEMTable(sem_mtx)
    res_a <- .defineNucleotideColors(
        sem_mtx_long = sem_mtx_long,
        motifSeq = "AAAAAAA"
    )
    # make sure we have a mseq and text_color column added
    expect_equal(colnames(res_a), c(
        "bp", "motif_pos",
        "sem_score", "mseq", "text_color"
    ))
    # all A nucleotides should be white
    expect_equal(
        res_a$text_color[res_a$bp == "A"],
        rep("white", 7)
    )
    # all other nucleotides should be "#d7dbdd" (grey)
    expect_equal(
        res_a$text_color[res_a$bp != "A"],
        rep("#d7dbdd", 21)
    )
    # mseq is A when bp is A
    expect_equal(
        res_a$mseq[res_a$bp == "A"],
        rep("A", 7)
    )
    # mseq is empty when bp is not A
    expect_equal(
        res_a$mseq[res_a$bp != "A"],
        rep("", 21)
    )
})


test_that("plotSEM returns a ggplot object", {
    plt_a <- plotSEM(SEMC,
        motif = "MA0099.2_HeLa",
        motifSeq = "AAAAAAA",
        highlight = 2
    )
    expect_s3_class(plt_a, "gg")
})
