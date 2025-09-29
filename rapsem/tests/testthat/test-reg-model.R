test_that("reg_model initializes correctly with simple model", {
  reg_lines <- c("med ~ treat + cov2",
                 "out ~ treat + med + cov1")
  model <- reg_model$new(reg_lines, intercept_name = "intercept")

  expect_s3_class(model, "reg_model")
  expect_equal(model$text, reg_lines)
  expect_named(model$med_preds, "med")
  expect_named(model$out_preds, "out")
  expect_equal(model$med_preds$med, c("intercept", "treat", "cov2"))
  expect_equal(model$out_preds$out, c("intercept", "treat", "med", "cov1"))
})

test_that("reg_model initializes correctly with split lines", {
  reg_lines <- c("med ~ treat",
                 "med ~ cov2",
                 "med ~ treat:cov2",
                 "out ~ treat",
                 "out ~ med",
                 "out ~ cov1")
  model <- reg_model$new(reg_lines, intercept_name = "intercept")

  expect_s3_class(model, "reg_model")
  expect_equal(model$text, reg_lines)
  expect_named(model$med_preds, "med")
  expect_named(model$out_preds, "out")
  expect_equal(model$med_preds$med, c("intercept", "treat", "cov2",
                                      "treat:cov2"))
  expect_equal(model$out_preds$out, c("intercept", "treat", "med", "cov1"))
})

test_that("reg_model identifies only outcomes correctly", {
  reg_lines <- c("out ~ treat + cov1")
  model <- reg_model$new(reg_lines, intercept_name = "intercept")

  expect_equal(names(model$med_preds), NULL)
  expect_named(model$out_preds, "out")
})

test_that("reg_model handles multiple mediators correctly", {
  reg_lines <- c(
    "med1 ~ treat + cov1",
    "med2 ~ treat + cov2",
    "out ~ treat + med1 + med2"
  )
  model <- reg_model$new(reg_lines, intercept_name = "intercept")

  expect_named(model$med_preds, c("med1", "med2"))
  expect_named(model$out_preds, "out")

  expect_equal(model$med_preds$med1, c("intercept", "treat", "cov1"))
  expect_equal(model$med_preds$med2, c("intercept", "treat", "cov2"))
  expect_equal(model$out_preds$out, c("intercept", "treat", "med1", "med2"))
})

test_that("reg_model handles multiple outcomes correctly", {
  reg_lines <- c(
    "out1 ~ treat + med",
    "out2 ~ treat + med + cov1",
    "med ~ treat"
  )
  model <- reg_model$new(reg_lines, intercept_name = "intercept")

  expect_named(model$med_preds, "med")
  expect_named(model$out_preds, c("out1", "out2"))

  expect_equal(model$med_preds$med, c("intercept", "treat"))
  expect_equal(model$out_preds$out1, c("intercept", "treat", "med"))
  expect_equal(model$out_preds$out2, c("intercept", "treat", "med", "cov1"))
})

test_that("reg_model handles zero intercepts correctly", {
  reg_lines <- c(
    "med ~ 0 + treat + cov1",
    "out ~ treat + med"
  )
  model <- reg_model$new(reg_lines, intercept_name = "intercept")

  expect_named(model$med_preds, "med")
  expect_named(model$out_preds, "out")
  expect_equal(model$med_preds$med, c("treat", "cov1"))
  expect_equal(model$out_preds$out, c("intercept", "treat", "med"))
})

test_that("reg_model adds main effects for interaction terms and warns", {
  reg_lines <- c(
    "out ~ treat + cov1 + cov1:cov2"
  )
  expect_warning(
    model <- reg_model$new(reg_lines, intercept_name = "intercept"),
    regexp = paste("Main effect cov2 added as predictor because interaction",
                   "cov1:cov2 was detected in outcome out")
  )
  expect_named(model$out_preds, "out")
  expect_true(all(c("intercept", "cov1", "cov2", "cov1:cov2")
                  %in% model$out_preds$out))
})
