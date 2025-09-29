test_that("var_name_config allows custom intercept name", {
  cfg <- var_name_config$new(treatment_names = "x",
                             all_var_names = c("x", "y"),
                             intercept_name = "const")
  expect_equal(cfg$intercept, "const")
})