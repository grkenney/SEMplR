test_that("SNPEffectMatrixCollection fails on invalid input", {
    expect_error(
        SNPEffectMatrixCollection(
            sems = getSEMs(SEMC)[[1]],
            semData = semData(SEMC)[1, ],
            semKey = ""
        ),
        regexp = "must provide a semKey if providing semData"
    )

    expect_error(
        SNPEffectMatrixCollection(
            sems = getSEMs(SEMC)[[1]],
            semData = semData(SEMC)[1, ],
            semKey = "foo"
        ),
        regexp = "semKey must be a column in semData"
    )
})


test_that("SNPEffectMatrixCollection remakes default SEMs collection", {
    semc_a <- SNPEffectMatrixCollection(
        sems = getSEMs(SEMC),
        semData = semData(SEMC),
        semKey = "transcription_factor"
    )
    expect_equal(semc_a, SEMC)

    m1 <- matrix(rnorm(16), nrow = 4)
    colnames(m1) <- c("A", "C", "G", "T")
    sem_list <- list(SNPEffectMatrix(sem = m1, baseline = 1, semId = "m1"))

    meta_data <- data.table::data.table(tf = "tf1", sem = "m1.sem")
    expect_message(
        SNPEffectMatrixCollection(
            sems = sem_list,
            semData = meta_data,
            semKey = "sem"
        ),
        regexp = "Removing .sem suffixes from semKey."
    )
})


test_that("SNPEffectMatrixCollection makes collection without .sem suffix", {
    meta_dt <- data.table::data.table(
        transcription_factor = "TFAP2B",
        SEM = "TFAP2B"
    )
    semc_a <- SNPEffectMatrixCollection(
        sems = getSEMs(SEMC)[[1]],
        semData = meta_dt,
        semKey = "SEM"
    )
    data.table::setkey(meta_dt, "SEM")
    # this test should match the first entry of the default data in the sems slot
    expect_equal(getSEMs(semc_a)[[1]], getSEMs(SEMC)[[1]])
    expect_equal(semData(semc_a), meta_dt)
    expect_equal(data.table::key(semData(semc_a)), "SEM")
})


test_that("sems pulls correct data types and lengths", {
    # pulling a single SEM by index
    sems_a <- getSEMs(SEMC)[[1]]
    expect_s4_class(sems_a, "SNPEffectMatrix")

    # pulling all SEMs
    sems_a <- getSEMs(SEMC)
    expect_type(sems_a, "list")
    expect_length(sems_a, 223)

    # pulling a single SEM by name
    sems_a <- getSEMs(SEMC, "ZSCAN4")
    expect_s4_class(sems_a, "SNPEffectMatrix")

    # pulling multiple SEMs by name
    sems_a <- getSEMs(SEMC, c("ZNF18", "ZSCAN4"))
    expect_type(sems_a, "list")
    expect_length(sems_a, 2)
})


test_that("show SNPEffectMatrixCollection collection", {
    expect_output(print(SEMC), "An object of class SNPEffectMatrixCollection")
    expect_output(print(SEMC), "SEMs\\(223\\): TFAP2B, ARNT")
    expect_output(print(SEMC), "semData\\(12\\): transcription_factor, ")
})
