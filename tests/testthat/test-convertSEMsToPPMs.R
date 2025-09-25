test_that("convertSEMsToPPMs errors on invalid input", {
    expect_error(convertSEMsToPPMs("foo"),
        regexp = "x must be of class SNPEffectMatrixCollection or"
    )
})


test_that("convertSEMsToPPMs handle different input types", {
    # SNPEffectMatrixCollection
    expect_no_condition(convertSEMsToPPMs(SEMC))

    # SNPEffectMatrix
    expect_no_condition(convertSEMsToPPMs(getSEMs(SEMC)[[1]]))

    # list of SNPEffectMatrices
    expect_no_condition(convertSEMsToPPMs(getSEMs(SEMC)[seq_len(2)]))
})
