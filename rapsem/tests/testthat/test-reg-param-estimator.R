source(test_path("helpers/mock_helpers.R"))

test_that("reg_param_estimator initializes correctly", {
  estimator <- reg_param_estimator$new(
    sample = mock_data,
    var_names = mock_var_names_config$new(),
    target_name = "out",
    factor_scores = mock_factor_scores$new(),
    predictor_names = c("intercept", "treat", "cov1", "cov2",
                        "cov1:cov2", "med", "cov2:med")
  )

  expect_true(!is.null(estimator$target))
  expect_equal(length(estimator$target), nrow(mock_data))
  expect_true(is.matrix(estimator$predictors))
  expect_equal(ncol(estimator$predictors), 7)
  expect_equal(colnames(estimator$predictors),
               c("intercept", "treat", "cov1", "cov2",
                 "cov1:cov2", "med", "cov2:med"))
  expect_true(is.matrix(estimator$mat_m1))
  expect_equal(dim(estimator$mat_m1), c(7, 7))
  expect_true(is.matrix(estimator$mat_m2))
  expect_equal(dim(estimator$mat_m2), c(7, 7))
  expect_true(is.numeric(estimator$vec_m1))
  expect_equal(length(estimator$vec_m1), 7)
  expect_true(is.numeric(estimator$vec_m2))
  expect_equal(length(estimator$vec_m2), 7)
})

test_that("reg_param_estimator performs regression and stores parameters", {
  estimator <- reg_param_estimator$new(
    sample = mock_data,
    var_names = mock_var_names_config$new(),
    target_name = "out",
    factor_scores = mock_factor_scores$new(),
    predictor_names = c("intercept", "treat", "cov1", "cov2",
                        "cov1:cov2", "med", "cov2:med")
  )

  expect_null(estimator$parameters)
  estimator$regress(FALSE, 0)
  expect_true(!is.null(estimator$parameters))
  expect_true(is.numeric(estimator$parameters))
  expect_equal(length(estimator$parameters), 7)
})

test_that("regression falls back to ginv on singular matrix", {
  estimator <- reg_param_estimator$new(
    sample = mock_data,
    var_names = mock_var_names_config$new(),
    target_name = "out",
    factor_scores = mock_factor_scores$new(),
    predictor_names = c("treat", "cov2", "med", "cov2_dup")
  )

  expect_warning(
    estimator$regress(FALSE, 0),
    "Using ginv due to singular matrix."
  )
  expect_true(!is.null(estimator$parameters))
  expect_true(is.numeric(estimator$parameters))
  expect_equal(length(estimator$parameters), 4)
})

test_that("reg_param_estimator calculates correct parameter estimates", {
  estimator <- reg_param_estimator$new(
    sample = mock_data,
    var_names = mock_var_names_config$new(),
    target_name = "out",
    factor_scores = mock_factor_scores$new(),
    predictor_names = c("intercept", "treat", "cov1", "med", "cov1:treat")
  )
  estimator$regress(FALSE, 0)
  expected_results <- get_regression_parameters()
  expect_equal(estimator$mat_m1, expected_results$mat_m1, tolerance = 1e-6)
  expect_equal(estimator$mat_m2, expected_results$mat_m2, tolerance = 1e-6)
  expect_equal(estimator$vec_m1, expected_results$vec_m1, tolerance = 1e-6)
  expect_equal(estimator$vec_m2, expected_results$vec_m2, tolerance = 1e-6)
  expect_equal(estimator$parameters, expected_results$theta, tolerance = 1e-6)
})

test_that("modify_moments applies correct  moment modification", {
  estimator <- reg_param_estimator$new(
    sample = mock_data,
    var_names = mock_var_names_config$new(),
    target_name = "out",
    factor_scores = mock_factor_scores$new(),
    predictor_names = c("intercept", "treat", "cov1", "med", "cov1:treat")
  )
  moments <- estimator$.__enclos_env__$private$modify_moments(
    estimator$mat_m1,
    estimator$vec_m1,
    estimator$mat_m2,
    estimator$vec_m2,
    estimator$target,
    estimator$var_target_error,
    modified_2smm = TRUE,
    tau = 2
  )
  r1 <- rbind(
    cbind(estimator$mat_m1, estimator$vec_m1),
    cbind(t(estimator$vec_m1), mean(estimator$target^2))
  )
  r2 <- rbind(
    cbind(-estimator$mat_m2, -estimator$vec_m2),
    cbind(t(estimator$vec_m2), estimator$var_target_error)
  )
  r1_inv <- solve(chol(r1))
  eig <- eigen(r1_inv %*% r2 %*% t(r1_inv), only.values = TRUE)$values
  lambda_max <- max(Re(eig))
  num_obs <- length(estimator$target)
  expected_factor <- if (lambda_max >= (1 + 1 / num_obs)) {
    1 - 2 / num_obs
  } else {
    1 / lambda_max - 1 / num_obs - 2 / num_obs
  }
  expected_mat <- estimator$mat_m1 + expected_factor * estimator$mat_m2
  expected_vec <- estimator$vec_m1 + expected_factor * estimator$vec_m2
  expect_equal(moments$mat, expected_mat, tolerance = 1e-8)
  expect_equal(moments$vec, expected_vec, tolerance = 1e-8)
})
