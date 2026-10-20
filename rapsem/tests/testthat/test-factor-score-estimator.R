source(test_path("helpers/mock_helpers.R"))

test_that("factor_score_estimator initializes and returns factor scores", {
  model <- mock_cfa_model$new()
  estimator <- factor_score_estimator$new(model, mock_data, FALSE)
  expect_s3_class(estimator, "factor_score_estimator")
  expect_true(!is.null(estimator$vals))
  expect_equal(nrow(estimator$vals), nrow(mock_data))
  expect_equal(ncol(estimator$vals), length(model$fact_names))
  expect_equal(sort(colnames(estimator$vals)), sort(model$fact_names))
  expect_true(!is.null(estimator$sigma))
  expect_equal(dim(estimator$sigma), c(3, 3))
  expect_true(isSymmetric(estimator$sigma))
  # default: Wall-Amemiya scores and Sigma_ee in the anchor metric, unscaled
  raw <- get_estimated_factor_scores(raw = TRUE)
  expect_false(estimator$standardized)
  expect_equal(estimator$vals, raw$factor_scores, tolerance = 1e-6)
  expect_equal(estimator$sigma, raw$sigma_ee, tolerance = 1e-6)
})

test_that("standardize = TRUE gives the latent-standardized metric", {
  estimator <- factor_score_estimator$new(mock_cfa_model$new(), mock_data,
                                          FALSE, standardize = TRUE)
  expected_results <- get_estimated_factor_scores()
  expect_true(estimator$standardized)
  expect_equal(estimator$vals, expected_results$factor_scores, tolerance = 1e-6)
  expect_equal(estimator$sigma, expected_results$sigma_ee, tolerance = 1e-6)
  latent_var <- apply(estimator$vals, 2, var) - diag(estimator$sigma)
  expect_equal(unname(colMeans(estimator$vals)), rep(0, 3), tolerance = 1e-10)
  expect_equal(unname(latent_var), rep(1, 3), tolerance = 1e-10)
})

test_that("estimator handles multiple missing values across indicators", {
  model <- mock_cfa_model$new()
  data_na <- mock_data
  data_na[sample(1:num_obs, 10), "i3"] <- NA
  data_na[sample(1:num_obs, 10), "i6"] <- NA
  estimator <- factor_score_estimator$new(model, data_na, FALSE)
  expect_equal(nrow(estimator$vals), num_obs)
})

test_that("estimator throws error when model does not converge", {
  model <- mock_cfa_model$new()
  bad_data <- data.frame(matrix(rnorm(100 * 9), ncol = 9))
  names(bad_data) <- paste0("i", 1:9)
  expect_error(factor_score_estimator$new(model, bad_data, FALSE),
               regexp = "CFA model did not converge")
})

test_that("pattern scores equal Wall-Amemiya scores on complete data", {
  raw <- get_estimated_factor_scores(raw = TRUE)
  private <- factor_score_estimator$new(mock_cfa_model$new(), mock_data,
                                        FALSE)$.__enclos_env__$private
  pattern <- private$pattern_scores(raw$z, rbind(raw$beta1hat, diag(3)),
                                    c(raw$beta0hat, rep(0, 3)), raw$psihat)
  expect_equal(unname(pattern$scores), unname(raw$factor_scores),
               tolerance = 1e-10)
  expect_equal(unname(pattern$sigma_rows[1, , ]), unname(raw$sigma_ee),
               tolerance = 1e-10)
})

test_that("missing indicator values give row-specific error covariances", {
  data_na <- mock_data
  data_na$i5[1:10] <- NA
  data_na$i2[11:15] <- NA
  estimator <- factor_score_estimator$new(mock_cfa_model$new(), data_na,
                                          FALSE)
  expect_true(all(estimator$rows))
  expect_equal(dim(estimator$sigma_rows), c(num_obs, 3, 3))
  # fewer observed indicators -> larger error variance of that factor
  expect_gt(estimator$sigma_rows[1, "med", "med"],
            estimator$sigma_rows[20, "med", "med"])
  expect_gt(estimator$sigma_rows[11, "cov1", "cov1"],
            estimator$sigma_rows[20, "cov1", "cov1"])
  # standardized: average error variance, row-specific values rescaled alike
  std <- factor_score_estimator$new(mock_cfa_model$new(), data_na, FALSE,
                                    standardize = TRUE)
  latent_var <- apply(std$vals, 2, var) - diag(std$sigma)
  expect_equal(unname(latent_var), rep(1, 3), tolerance = 1e-10)
  expect_equal(std$sigma_rows[1, , ] / estimator$sigma_rows[1, , ],
               std$sigma / estimator$sigma, tolerance = 1e-10)
})

test_that("rows without factor indicators or structural data are dropped", {
  data_na <- mock_data
  data_na[1:3, c("i4", "i5", "i6")] <- NA
  data_na$cov2[4] <- NA
  estimator <- factor_score_estimator$new(mock_cfa_model$new(), data_na,
                                          FALSE, keep_vars = "cov2")
  expect_equal(which(!estimator$rows), 1:4)
  expect_equal(nrow(estimator$vals), num_obs - 4)
  expect_equal(dim(estimator$sigma_rows)[1], num_obs - 4)
})

test_that("reliability is McDonald's omega over all indicators", {
  estimator <- factor_score_estimator$new(mock_cfa_model$new(), mock_data,
                                          FALSE)
  fit <- lavaan::cfa(mock_cfa_model$new()$lavaan_syntax, data = mock_data,
                     meanstructure = TRUE, estimator = "ML", se = "none")
  pe <- lavaan::parameterEstimates(fit)
  omega <- sapply(latent_names, function(f) {
    rows <- pe$lhs == f & pe$op == "=~"
    theta <- pe$est[pe$op == "~~" & pe$lhs %in% pe$rhs[rows] &
                      pe$lhs == pe$rhs]
    psi <- pe$est[pe$op == "~~" & pe$lhs == f & pe$rhs == f]
    sum(pe$est[rows])^2 * psi / (sum(pe$est[rows])^2 * psi + sum(theta))
  })
  expect_equal(estimator$reliability, omega, tolerance = 1e-6)
})

test_that("blockwise estimation equals the orthogonal joint fit", {
  blockwise <- factor_score_estimator$new(mock_cfa_model$new(), mock_data,
                                          FALSE, blockwise = TRUE)
  orthogonal <- factor_score_estimator$new(mock_cfa_model$new(), mock_data,
                                           TRUE)
  expect_equal(blockwise$vals, orthogonal$vals, tolerance = 1e-4)
  expect_equal(blockwise$sigma, orthogonal$sigma, tolerance = 1e-4)
  expect_equal(blockwise$reliability, orthogonal$reliability,
               tolerance = 1e-4)
})

test_that("blockwise estimation keeps within-block residual covariances", {
  model <- mock_cfa_model$new()
  model$lavaan_syntax <- paste(model$lavaan_syntax, "i5 ~~ i6", sep = "\n")
  estimator <- factor_score_estimator$new(model, mock_data, FALSE,
                                          blockwise = TRUE)
  expect_equal(dim(estimator$vals), c(num_obs, 3))
})

test_that("blockwise estimation rejects cross-block residual covariances", {
  model <- mock_cfa_model$new()
  model$lavaan_syntax <- paste(model$lavaan_syntax, "i2 ~~ i5", sep = "\n")
  expect_error(factor_score_estimator$new(model, mock_data, FALSE,
                                          blockwise = TRUE),
               "different latent variables")
})

test_that("blockwise estimation requires three indicators per factor", {
  model <- mock_cfa_model$new()
  model$lavaan_syntax <- "cov1 =~ i1 + i2\nmed =~ i4 + i5 + i6\nout =~ i7 + i8 + i9"
  model$free_ind$cov1 <- "i2"
  expect_error(factor_score_estimator$new(model, mock_data, FALSE,
                                          blockwise = TRUE),
               "fewer than 3 indicators")
})

heywood_data <- function() {
  # one almost error-free indicator: sampling error gives a negative, but not
  # significantly negative, residual variance (-0.62, SE 3.1)
  set.seed(3)
  n <- 60
  f <- rnorm(n)
  data.frame(i1 = f + rnorm(n, sd = 0.05), i2 = 0.5 * f + rnorm(n),
             i3 = 0.5 * f + rnorm(n))
}

improper_data <- function(n = 1000) {
  # misspecified one-factor model (correlated errors of i2 and i3): the
  # residual variance of i1 is significantly negative (z = -2.4)
  set.seed(1)
  f <- rnorm(n)
  e <- rnorm(n)
  data.frame(i1 = 0.8 * f + rnorm(n, sd = 0.6), i2 = 0.7 * f + e,
             i3 = 0.7 * f - 0.35 * e + rnorm(n, sd = 0.8))
}

test_that("clearly improper CFA solutions (Heywood cases) are rejected", {
  model <- cfa_model$new("f =~ i1 + i2 + i3")
  expect_error(factor_score_estimator$new(model, improper_data(), FALSE),
               "variance of i1 significantly negative")
  expect_error(factor_score_estimator$new(model, improper_data(), FALSE,
                                          blockwise = TRUE),
               "variance of i1 significantly negative")
})

test_that("non-significant negative residual variances are fixed at 0", {
  model <- cfa_model$new("f =~ i1 + i2 + i3")
  data <- heywood_data()
  fs <- factor_score_estimator$new(model, data, FALSE)
  # the error-free indicator determines the score; no estimation error
  expect_equal(cor(fs$vals[, 1], data$i1), 1)
  expect_equal(unname(fs$sigma[1, 1]), 0)
  expect_equal(factor_score_estimator$new(model, data, FALSE,
                                          blockwise = TRUE)$vals, fs$vals)
  # pattern-wise scores reach the same limit
  data$i2[1:5] <- NA
  fs_na <- factor_score_estimator$new(model, data, FALSE)
  expect_equal(fs_na$vals[6:60, 1], fs$vals[6:60, 1], tolerance = 1e-6)
  expect_lt(fs_na$sigma[1, 1], 1e-6)
})

test_that("bootstrap replicates with improper solutions are skipped", {
  good <- mock_data[, names(mock_data) != "intercept"]
  model <- "cov1 =~ i1 + i2 + i3
            med =~ i4 + i5 + i6
            out =~ i7 + i8 + i9
            med ~ treat + cov1 + treat:cov1
            out ~ treat + med + cov1"
  data_info <- data_loader$new(good, "treat")
  models <- model_loader$new(model, data_info$name_config)
  # replace the second bootstrap sample by one whose cov1 block is improper
  samples <- data_info$bootstrap(3)
  bad <- samples[[2]][rep(seq_len(nrow(samples[[2]])), length.out = 1000), ]
  bad[, c("i1", "i2", "i3")] <- improper_data()
  samples[[2]] <- bad
  expect_warning(
    res <- run_bootstrap_reps(3, samples, models, data_info, TRUE, 5, 1e-4,
                              FALSE, blockwise = TRUE),
    "1 of 3 bootstrap replicates failed"
  )
  expect_length(res$results_gestimation, 2)
})
