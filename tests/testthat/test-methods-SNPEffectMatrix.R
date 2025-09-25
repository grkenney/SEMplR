test_that("getSEM gets the SEM", {
    s <- getSEMs(SEMC, "MA0099.2_HeLa")
    sm_a <- getSEM(s)
    # spot check first two positions
    sm_e <- matrix(
        c(
            -0.985316, -1.2757100, -1.3195100, 0.0000000,
            -0.697860, -0.7446700, 0.0162566, -0.4440620
        ),
        ncol = 4, byrow = TRUE
    ) |>
        data.table::as.data.table() |>
        stats::setNames(c("A", "C", "G", "T"))
    expect_equal(sm_a[seq_len(2), ], sm_e, tolerance = 1e-5)
})


test_that("getBaseline gets the SEM baseline", {
    s <- getSEMs(SEMC, "MA0099.2_HeLa")
    expect_equal(getBaseline(s), -0.933424,
        tolerance = 1e-5
    )
})


test_that("getSEMId gets the SEM Id", {
    s <- getSEMs(SEMC, "MA0099.2_HeLa")
    expect_equal(getSEMId(s), "MA0099.2_HeLa")
})


test_that("show SNPEffectMatrix prints correctly", {
    m <- matrix(c(0, 0, 0, 0),
        ncol = 4,
        dimnames = list(c(1), c("A", "C", "G", "T"))
    )
    sem <- SNPEffectMatrix(
        sem = m,
        baseline = 1, semId = "SEM_ID"
    )
    expect_output(print(sem), "An object of class SNPEffectMatrix")
    expect_output(print(sem), "semId:  SEM_ID")
    expect_output(print(sem), "baseline:  1")
    expect_output(print(sem), "1:     0     0     0     0")
})
