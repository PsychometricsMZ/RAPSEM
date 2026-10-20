source(test_path("helpers/mock_helpers.R"))

test_that("summarize_results groups replicates by outcome name", {
  results <- list(out_rep1 = c(a = 1), med2_rep1 = c(a = 10),
                  out_rep2 = c(a = 3), med2_rep3 = c(a = 30))
  summary <- summarize_results(results, c("out", "med2"))
  expect_equal(summary$out["a", "Estimate"], 2)
  expect_equal(summary$med2["a", "Estimate"], 20)
})

test_that("the lavaan comparison reports the treated group", {
  model <- "med =~ i4 + i5 + i6
            out =~ i7 + i8 + i9
            cov1 =~ i1 + i2 + i3
            med ~ treat + cov1 + treat:cov1
            out ~ treat + med + cov1"
  vars <- rapsem:::data_loader$new(mock_data[, names(mock_data) != "intercept"],
                                   "treat")$name_config
  models <- model_loader$new(model, vars)
  estimate <- function(data) {
    fit <- lavaan::sem(models$group_lavaan, data = data, group = "treat",
                       group.equal = c("loadings", "intercepts"))
    summarize_lavaan_results(list(treat = fit), models, c(treat = "1"))$out
  }
  treated_first <- estimate(mock_data[order(-mock_data$treat), ])
  control_first <- estimate(mock_data[order(mock_data$treat), ])
  expect_equal(treated_first$Estimate, control_first$Estimate,
               tolerance = 1e-4)
})

test_that("med_est runs with blockwise and joint measurement models", {
  model <- "cov1 =~ i1 + i2 + i3
            med =~ i4 + i5 + i6
            out =~ i7 + i8 + i9
            med ~ treat + cov1 + treat:cov1
            out ~ treat + med + cov1"
  data <- mock_data[, names(mock_data) != "intercept"]
  set.seed(1)
  blockwise <- med_est(model, data, "treat", nboot = 3)
  set.seed(1)
  orthogonal <- med_est(model, data, "treat", nboot = 3, blockwise = FALSE,
                        orthogonal = TRUE)
  expect_equal(blockwise$summary$gestimation$out["med", "Estimate"],
               orthogonal$summary$gestimation$out["med", "Estimate"],
               tolerance = 1e-3)
})

test_that("the standardized lavaan comparison reports std.all", {
  model <- "med =~ i4 + i5 + i6
            out =~ i7 + i8 + i9
            cov1 =~ i1 + i2 + i3
            med ~ treat + cov1 + treat:cov1
            out ~ treat + med + cov1"
  vars <- rapsem:::data_loader$new(mock_data[, names(mock_data) != "intercept"],
                                   "treat")$name_config
  models <- model_loader$new(model, vars)
  fit <- lavaan::sem(models$group_lavaan, data = mock_data, group = "treat",
                     group.equal = c("loadings", "intercepts"))
  std <- summarize_lavaan_results(list(treat = fit), models, c(treat = "1"),
                                  standardize = TRUE)$out
  ref <- lavaan::standardizedSolution(fit, type = "std.all")
  group <- match("1", lavaan::lavInspect(fit, "group.label"))
  ref <- ref[ref$lhs == "out" & ref$op == "~" & ref$group == group, ]
  expect_equal(std$Estimate, ref$est.std)
})
