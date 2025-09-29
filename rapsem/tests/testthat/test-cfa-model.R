test_that("Indicators are correctly sorted", {
  syntax <- "
    med =~ i1 + i2 + i3
    out =~ i4 + i5 +i6
  "
  model <- cfa_model$new(syntax)
  testthat::expect_equal(model$fixed_ind$med, c("i1"))
  testthat::expect_equal(model$fixed_ind$out, c("i4"))
  testthat::expect_equal(model$free_ind$med, c("i2", "i3"))
  testthat::expect_equal(model$free_ind$out, c("i5", "i6"))
})
