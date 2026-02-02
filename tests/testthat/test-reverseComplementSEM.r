test_that("reverseComplementSEM generates correct output on JUN motif", {
  jun_sem <- getSEMs(SEMC, "JUN")
  rcs_a <- reverseComplementSEM(jun_sem)
    
  expect_s4_class(rcs_a, "SNPEffectMatrix")
  expect_equal(getBaseline(rcs_a), -0.933424)
  expect_equal(getSEM(rcs_a)$`T`, rev(getSEM(jun_sem)$A))
  expect_equal(getSEM(rcs_a)$A, rev(getSEM(jun_sem)$`T`))
  expect_equal(getSEM(rcs_a)$C, rev(getSEM(jun_sem)$G))
  expect_equal(getSEM(rcs_a)$G, rev(getSEM(jun_sem)$C))
})


test_that("reverseComplementSEM fails on invalid input", {
  expect_error(reverseComplementSEM(1), 
               "sem must be an object of class SNPEffectMatrix")
})
