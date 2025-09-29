source(test_path("helpers/mock_helpers.R"))

test_that("factor_score_estimator initializes and returns factor scores", {
  model <- mock_cfa_model$new()
  estimator <- factor_score_estimator$new(model, mock_data)
  expected_results <- get_estimated_factor_scores()
  expect_s3_class(estimator, "factor_score_estimator")
  expect_true(!is.null(estimator$vals))
  expect_equal(nrow(estimator$vals), nrow(mock_data))
  expect_equal(ncol(estimator$vals), length(model$fact_names))
  expect_equal(sort(colnames(estimator$vals)), sort(model$fact_names))
  expect_equal(estimator$vals, expected_results$factor_scores, tolerance = 1e-6)
  expect_true(!is.null(estimator$sigma))
  expect_equal(dim(estimator$sigma), c(3, 3))
  expect_true(isSymmetric(estimator$sigma))
  expect_equal(estimator$sigma, expected_results$sigma_ee, tolerance = 1e-6)
})

test_that("estimator handles multiple missing values across indicators", {
  model <- mock_cfa_model$new()
  data_na <- mock_data
  data_na[sample(1:num_obs, 10), "i3"] <- NA
  data_na[sample(1:num_obs, 10), "i6"] <- NA
  estimator <- factor_score_estimator$new(model, data_na)
  expect_equal(nrow(estimator$vals), num_obs)
})

test_that("estimator throws error when model does not converge", {
  model <- mock_cfa_model$new()
  bad_data <- data.frame(matrix(rnorm(100 * 9), ncol = 9))
  names(bad_data) <- paste0("i", 1:9)
  expect_error(factor_score_estimator$new(model, bad_data),
               regexp = "CFA model did not converge")
})
