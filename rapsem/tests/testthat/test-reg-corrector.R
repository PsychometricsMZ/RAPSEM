source(test_path("helpers/mock_helpers.R"))

test_that("reg_corrector calculates correct parameter estimates", {
  med_predictors <- list(med = c("intercept", "treat", "cov1", "cov1:treat"))
  out_predictors <- list(out = pred_names)
  corrector <- reg_corrector$new(
    reg_est = mock_reg_param_estimator$new(),
    reg_model = mock_reg_model$new(med_preds = med_predictors,
                                   out_preds = out_predictors),
    sample = mock_data,
    var_names = mock_var_names_config$new(),
    factor_scores = mock_factor_scores$new(),
    modified_2smm = FALSE,
    tau = 0
  )
  corr_est <- corrector$corrected_reg_estimator
  expected_results <- get_corrected_reg_parms(1e-4)
  expect_equal(corr_est$mat_m1, expected_results$mat_m1, tolerance = 1e-6)
  expect_equal(corr_est$mat_m2, expected_results$mat_m2, tolerance = 1e-6)
  expect_equal(corr_est$vec_m1, expected_results$vec_m1, tolerance = 1e-6)
  expect_equal(corr_est$vec_m2, expected_results$vec_m2, tolerance = 1e-6)
  corrector$g_estimation(FALSE, 0, 1e-4)
  expect_equal(corr_est$parameters, expected_results$theta, tolerance = 1e-6)
})
