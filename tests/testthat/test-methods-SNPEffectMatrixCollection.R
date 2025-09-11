test_that("SNPEffectMatrixCollection fails on invalid input", {
  expect_error(SNPEffectMatrixCollection(sems = sems(sc)[[1]],
                                         semData = semData(sc)[1, ],
                                         semKey = ""),
               regexp = "must provide a semKey if providing semData")

  expect_error(SNPEffectMatrixCollection(sems = sems(sc)[[1]],
                                         semData = semData(sc)[1, ],
                                         semKey = "foo"),
               regexp = "semKey must be a column in semData")
})


test_that("SNPEffectMatrixCollection remakes default SEMs collection", {
  semc_a <- SNPEffectMatrixCollection(sems = sems(sc),
                                      semData = semData(sc),
                                      semKey = "SEM")
  expect_equal(semc_a, sc)
  expect_message(SNPEffectMatrixCollection(sems = sems(sc),
                                           semData = semData(sc),
                                           semKey = "SEM"),
                 regexp = "Removing .sem suffixes from semKey.")
})


test_that("SNPEffectMatrixCollection makes collection without .sem suffix", {
  meta_dt <- data.table::data.table(transcription_factor = "TFAP2B",
                                    SEM = "AP2B_HUMAN.SK-N-SH")
  semc_a <- SNPEffectMatrixCollection(sems = sems(sc)[[1]],
                                      semData = meta_dt,
                                      semKey = "SEM")
  data.table::setkey(meta_dt, "SEM")
  # this test should match the first entry of the default data in the sems slot
  expect_equal(sems(semc_a)[[1]], sems(sc)[[1]])
  expect_equal(semData(semc_a), meta_dt)
  expect_equal(data.table::key(semData(semc_a)), "SEM")
})


test_that("sems pulls correct data types and lengths", {
  # pulling a single SEM by index
  sems_a <- sems(sc)[[1]]
  expect_s4_class(sems_a, "SNPEffectMatrix")
  
  # pulling all SEMs
  sems_a <- sems(sc)
  expect_type(sems_a, "list")
  expect_length(sems_a, 223)
  
  # pulling a single SEM by name
  sems_a <- sems(sc, "ZSCAN4_secondary")
  expect_s4_class(sems_a, "SNPEffectMatrix")
  
  # pulling multiple SEMs by name
  sems_a <- sems(sc, c("ZNF18_HUMAN.HEK293", "ZSCAN4_secondary"))
  expect_type(sems_a, "list")
  expect_length(sems_a, 2)
})


test_that("show SNPEffectMatrixCollection collection", {
  expect_output(print(sc), "An object of class SNPEffectMatrixCollection")
  expect_output(print(sc), "sems\\(223\\): AP2B_HUMAN.SK-N-SH, ")
  expect_output(print(sc), "semData\\(13\\): transcription_factor, ")
})

