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
  expect_equal(estimator$sensitivity_param,
               expected_results$sensitivity_param, tolerance = 1e-6)
})

# factor of the modified 2SMM estimator for a given estimator object
modify_factor <- function(estimator, tau) {
  estimator$.__enclos_env__$private$get_factor_to_modify_moments(
    estimator$mat_m1, estimator$vec_m1, estimator$mat_m2, estimator$vec_m2,
    estimator$target, estimator$var_target_error,
    modified_2smm = TRUE, tau = tau
  )
}

# reference implementation of Wall & Amemiya (2000, eqs. 23-24), written
# independently of the package code (Cholesky whitening instead of R1^(-1/2))
wall_amemiya_factor <- function(estimator, tau) {
  r1 <- rbind(
    cbind(estimator$mat_m1, estimator$vec_m1),
    cbind(t(estimator$vec_m1), mean(estimator$target^2))
  )
  r2 <- rbind(
    cbind(-estimator$mat_m2, -estimator$vec_m2),
    cbind(-t(estimator$vec_m2), estimator$var_target_error)
  )
  u_inv <- solve(chol(r1))
  gamma_hat <- max(eigen(t(u_inv) %*% r2 %*% u_inv, symmetric = TRUE,
                         only.values = TRUE)$values)
  num_obs <- length(estimator$target)
  if (1 / gamma_hat >= 1 + 1 / num_obs) {
    1 - tau / num_obs
  } else {
    1 / gamma_hat - 1 / num_obs - tau / num_obs
  }
}

new_out_estimator <- function(data = mock_data,
                              factor_scores = mock_factor_scores$new()) {
  reg_param_estimator$new(
    sample = data,
    var_names = mock_var_names_config$new(),
    target_name = "out",
    factor_scores = factor_scores,
    predictor_names = c("intercept", "treat", "cov1", "med", "cov1:treat")
  )
}

scaled_error_scores <- function(scale) {
  sigma <- scale * diag(3)
  dimnames(sigma) <- list(latent_names, latent_names)
  R6::R6Class("factor_score_estimator", public = list(
    vals = latent_vars,
    sigma = sigma,
    names = latent_names
  ))$new()
}

test_that("modify_moments follows Wall & Amemiya (2000, eq. 24)", {
  estimator <- new_out_estimator()
  expect_equal(modify_factor(estimator, tau = 2),
               wall_amemiya_factor(estimator, tau = 2), tolerance = 1e-8)
  expect_lte(modify_factor(estimator, tau = 2), 1)
})

test_that("modify_moments uses the full correction when it is safe", {
  estimator <- new_out_estimator(factor_scores = scaled_error_scores(0.01))
  num_obs <- length(estimator$target)
  expect_equal(modify_factor(estimator, tau = 5), 1 - 5 / num_obs)
})

test_that("modify_moments shrinks the correction near singularity", {
  estimator <- new_out_estimator(factor_scores = scaled_error_scores(0.9))
  num_obs <- length(estimator$target)
  factor <- modify_factor(estimator, tau = 5)
  expect_lt(factor, 1 - 5 / num_obs)
  expect_equal(factor, wall_amemiya_factor(estimator, tau = 5),
               tolerance = 1e-8)
})

test_that("modify_moments is 1 - tau/n without measurement error", {
  estimator <- new_out_estimator(factor_scores = scaled_error_scores(0))
  num_obs <- length(estimator$target)
  expect_equal(modify_factor(estimator, tau = 5), 1 - 5 / num_obs)
})

test_that("modify_moments does not depend on the treatment coding", {
  data_pm <- mock_data
  data_pm$treat <- 2 * data_pm$treat - 1
  for (scale in c(0.01, 0.9)) {
    scores <- scaled_error_scores(scale)
    expect_equal(modify_factor(new_out_estimator(mock_data, scores), 5),
                 modify_factor(new_out_estimator(data_pm, scores), 5),
                 tolerance = 1e-8)
  }
})

test_that("modify_moments uses R1^-1 R2 for the g-estimation system", {
  # identified RAPSEM outcome equation: no treatment interaction in the
  # outcome, so the mediator instrument is not collinear with other rows
  out_preds <- c("intercept", "treat", "cov1", "med")
  out_est <- reg_param_estimator$new(
    sample = mock_data,
    var_names = mock_var_names_config$new(),
    target_name = "out",
    factor_scores = mock_factor_scores$new(),
    predictor_names = out_preds
  )
  corrector <- reg_corrector$new(
    reg_est = out_est,
    reg_model = mock_reg_model$new(
      med_preds = list(med = c("intercept", "treat", "cov1", "cov1:treat")),
      out_preds = list(out = out_preds)),
    sample = mock_data,
    var_names = mock_var_names_config$new(),
    factor_scores = mock_factor_scores$new(),
    modified_2smm = FALSE,
    tau = 0
  )
  g_est <- corrector$corrected_reg_estimator
  r1 <- rbind(cbind(g_est$mat_m1, g_est$vec_m1),
              cbind(t(g_est$vec_m1), mean(g_est$target^2)))
  r2 <- rbind(cbind(-g_est$mat_m2, -g_est$vec_m2),
              cbind(-t(g_est$vec_m2), g_est$var_target_error))
  expect_gt(max(abs(r1 - t(r1))), 1e-3)
  gamma_hat <- max(Re(eigen(solve(r1) %*% r2, only.values = TRUE)$values))
  num_obs <- length(g_est$target)
  expected <- if (1 / gamma_hat >= 1 + 1 / num_obs) {
    1 - 5 / num_obs
  } else {
    1 / gamma_hat - 1 / num_obs - 5 / num_obs
  }
  expect_equal(modify_factor(g_est, tau = 5), expected, tolerance = 1e-8)
})

test_that("interaction columns are the product of all their parts", {
  estimator <- reg_param_estimator$new(
    sample = mock_data,
    var_names = mock_var_names_config$new(),
    target_name = "out",
    factor_scores = mock_factor_scores$new(),
    predictor_names = c("intercept", "treat", "cov2", "cov1",
                        "treat:cov2:cov1")
  )
  x <- estimator$predictors
  expect_equal(x[, "treat:cov2:cov1"], x[, "treat"] * x[, "cov2"] * x[, "cov1"])
})

test_that("observed variables are centered, and scaled only if standardized", {
  private <- new_out_estimator()$.__enclos_env__$private
  vals <- c(1, 2, 3, 10)
  expect_equal(private$center_if_necessary("z", vals, "treat", "intercept"),
               vals - mean(vals))
  expect_equal(as.numeric(private$center_if_necessary(
    "z", vals, "treat", "intercept", "cov1", standardize = TRUE)),
    as.numeric(scale(vals)))
  # factor scores: only centered (standardized by factor_score_estimator)
  expect_equal(private$center_if_necessary(
    "cov1", vals, "treat", "intercept", "cov1", standardize = TRUE),
    vals - mean(vals))
  expect_equal(private$center_if_necessary("treat", c(0, 1, 1), "treat",
                                           "intercept", standardize = TRUE),
               c(0, 1, 1))
})

test_that("standardized estimates are the rescaled anchor-metric estimates", {
  preds <- c("intercept", "treat", "cov1", "med", "cov1:treat")
  raw <- factor_score_estimator$new(mock_cfa_model$new(), mock_data, FALSE)
  std <- factor_score_estimator$new(mock_cfa_model$new(), mock_data, FALSE,
                                    standardize = TRUE)
  latent_sd <- sqrt(apply(raw$vals, 2, var) - diag(raw$sigma))
  # coefficient_std = coefficient_raw * sd(predictor) / sd(outcome)
  scale <- c(intercept = 1, treat = 1, cov1 = latent_sd[["cov1"]],
             med = latent_sd[["med"]], `cov1:treat` = latent_sd[["cov1"]]) /
    latent_sd[["out"]]
  for (modified in c(FALSE, TRUE)) {
    est <- lapply(list(raw = raw, std = std), function(fs) {
      reg_param_estimator$new(mock_data, mock_var_names_config$new(), "out",
                              fs, preds)$regress(modified, 5)[preds]
    })
    expect_equal(est$std, est$raw * scale[preds], tolerance = 1e-8)
  }
})
