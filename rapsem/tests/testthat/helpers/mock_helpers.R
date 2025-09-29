mock_var_names_config <- R6::R6Class("var_names_config", public = list(
  all = c(paste("i", 1:9, sep = ""), "out", "treat", "cov1", "cov2", "med"),
  intercept = "intercept",
  treatment = c("treat")
))

mock_cfa_model <- R6::R6Class("cfa_model", public = list(
  fact_names = c("cov1", "med", "out"),
  free_ind = list(
    cov1 = c("i2", "i3"),
    med = c("i5", "i6"),
    out = c("i8", "i9")
  ),
  fixed_ind = list(
    cov1 = c("i1"),
    med = c("i4"),
    out = c("i7")
  ),
  lavaan_syntax = "cov1 =~ i1 + i2 + i3
                med =~ i4 + i5 + i6
                out =~ i7 + i8 + i9",
  txt = NULL,
  initialize = function(txt = NULL) {
    self$txt <- txt
  }
))

mock_reg_model <- R6::R6Class("reg_model", public = list(
  txt = NULL,
  intercept = NULL,
  med_preds = NULL,
  out_preds = NULL,
  initialize = function(txt = NULL, intercept = NULL,
                        med_preds = NULL, out_preds = NULL) {
    self$txt <- txt
    self$intercept <- intercept
    self$med_preds <- med_preds
    self$out_preds <- out_preds
  }
))

set.seed(123)
num_obs <- 100

treat <- sample(0:1, num_obs, replace = TRUE)
cov1 <- rnorm(num_obs)
cov2 <- rnorm(num_obs)
med  <- (0.3 * treat + 0.2 * cov1 + 0.4 * cov1 * treat + 0.2 * cov2
         + rnorm(num_obs, sd = 0.5))
out  <- (0.2 * treat + 0.2 * cov1 + 0.2 * cov2 + 0.6 * med
         + rnorm(num_obs, sd = 0.6))

mock_data <- data.frame(
  i1 = cov1 + rnorm(num_obs, sd = 0.2),
  i2 = cov1 + rnorm(num_obs, sd = 0.2),
  i3 = cov1 + rnorm(num_obs, sd = 0.2),

  i4 = med + rnorm(num_obs, sd = 0.2),
  i5 = med + rnorm(num_obs, sd = 0.2),
  i6 = med + rnorm(num_obs, sd = 0.2),

  i7 = out + rnorm(num_obs, sd = 0.2),
  i8 = out + rnorm(num_obs, sd = 0.2),
  i9 = out + rnorm(num_obs, sd = 0.2),

  cov2 = cov2,
  cov2_dup = cov2,
  treat = treat,
  intercept = rep(1, num_obs)
)

get_estimated_factor_scores <- function() {
  data_transposed <- mock_data[, c(2, 3, 5, 6, 8, 9, 1, 4, 7)]
  model_cfa <- mock_cfa_model$new()$lavaan_syntax
  cfa <- lavaan::cfa(model_cfa, data = data_transposed,
                     meanstructure = TRUE, estimator = "MLR",
                     missing = "ML", se = "none", std.lv = TRUE)
  beta0hat <- matrix(lavaan::coef(cfa)[paste0("i", c(2, 3, 5, 6, 8, 9),
                                              "~1")], 6, 1)
  beta1hat <- matrix(0, 6, 3)
  beta1hat[1:2, 1] <- lavaan::coef(cfa)[paste0("cov1=~i", c(2, 3))]
  beta1hat[3:4, 2] <- lavaan::coef(cfa)[paste0("med=~i", c(5, 6))]
  beta1hat[5:6, 3] <- lavaan::coef(cfa)[paste0("out=~i", c(8, 9))]
  psihat <- matrix(0, 9, 9)
  diag(psihat) <- lavaan::coef(cfa)[paste0("i", c(2, 3, 5, 6, 8, 9, 1, 4, 7),
                                           "~~i", c(2, 3, 5, 6, 8, 9, 1, 4, 7))]
  b1_mod1 <- cbind(matrix(c(0), 3, 6), diag(3))
  b1_mod2 <- rbind(diag(6), (-1) * t(beta1hat))

  gamma <- (b1_mod1 %*% psihat %*% b1_mod2
            %*% solve(t(b1_mod2) %*% psihat %*% (b1_mod2)))
  sigma_ee <- (cbind(-gamma, (diag(3) + gamma %*% beta1hat))
               %*% psihat %*%  t(b1_mod1))

  factor_scores <- matrix(c(0), num_obs, 3)
  for (i in 1:num_obs){
    centered <- t(as.matrix(data_transposed[i, ])) - rbind(beta0hat,
                                                           matrix(0, 3, 1))
    factor_scores[i, ] <- (cbind(-gamma, (diag(3) + gamma %*% beta1hat))
                           %*% centered)
  }
  colnames(factor_scores) <- latent_names
  colnames(sigma_ee) <- rownames(sigma_ee) <- latent_names
  return(list(factor_scores = factor_scores, sigma_ee = sigma_ee))
}

latent_vars <- cbind(cov1 = cov1, med = med, out = out)
latent_cov <- cov(latent_vars)
latent_names <- colnames(latent_vars)
colnames(latent_cov) <- rownames(latent_cov) <- latent_names
mock_factor_scores <- R6::R6Class("factor_score_estimator", public = list(
  vals = latent_vars,
  sigma = latent_cov,
  names = latent_names
))

r0 <- mock_data[, "treat"]
x1 <- scale(latent_vars[, "cov1"])[, 1]
m  <- scale(latent_vars[, "med"])[, 1]
y  <- scale(latent_vars[, "out"])[, 1]
sigma <- latent_cov
preds <- cbind(
  intercept = mock_data[, "intercept"],
  treat = r0,
  cov1 = x1,
  med = m,
  `cov1:treat` = x1 * r0
)
pred_names <- colnames(preds)
target <- matrix(y, ncol = 1)
colnames(target) <- "out"

get_regression_parameters <- function() {
  bi11 <- 1
  bi12 <- mean(r0)
  bi13 <- mean(x1)
  bi14 <- mean(m)
  bi15 <- mean(x1 * r0)

  bi22 <- mean(r0^2)
  bi23 <- mean(r0 * x1)
  bi24 <- mean(r0 * m)
  bi25 <- mean(r0 * x1 * r0)

  bi33 <- mean(x1^2)
  bi34 <- mean(x1 * m)
  bi35 <- mean(x1^2 * r0)

  bi44 <- mean(m^2)
  bi45 <- mean(m * x1 * r0)

  bi55 <- mean(x1^2 * r0^2)

  mat_out_m1 <- matrix(c(
    bi11, bi12, bi13, bi14, bi15,
    bi12, bi22, bi23, bi24, bi25,
    bi13, bi23, bi33, bi34, bi35,
    bi14, bi24, bi34, bi44, bi45,
    bi15, bi25, bi35, bi45, bi55
  ), 5, 5, byrow = TRUE)
  colnames(mat_out_m1) <- rownames(mat_out_m1) <- pred_names

  ci33 <- mean(- sigma[1, 1])
  ci34 <- mean(- sigma[1, 2])
  ci35 <- mean(- sigma[1, 1] * r0)

  ci44 <- mean(- sigma[2, 2])
  ci45 <- mean(- sigma[1, 2] * r0)

  ci55 <- mean(- sigma[1, 1] * r0^2)

  mat_out_m2 <- matrix(c(
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, ci33, ci34, ci35,
    0, 0, ci34, ci44, ci45,
    0, 0, ci35, ci45, ci55
  ), 5, 5, byrow = TRUE)
  colnames(mat_out_m2) <- rownames(mat_out_m2) <- pred_names

  mat_out <- mat_out_m1 + mat_out_m2

  bi1 <- mean(y)
  bi2 <- mean(y * r0)
  bi3 <- mean(y * x1)
  bi4 <- mean(y * m)
  bi5 <- mean(y * x1 * r0)

  vec_out_m1 <- c(bi1, bi2, bi3, bi4, bi5)
  names(vec_out_m1) <- pred_names

  ci3 <- mean(- sigma[1, 3])
  ci4 <- mean(- sigma[2, 3])
  ci5 <- mean(- sigma[1, 3] * r0)

  vec_out_m2 <- c(0, 0, ci3, ci4, ci5)
  names(vec_out_m2) <- pred_names

  vec_out <- vec_out_m1 + vec_out_m2

  theta <- setNames(as.numeric(solve(mat_out, vec_out)), pred_names)

  list(theta = theta, mat_m1 = mat_out_m1, mat_m2 = mat_out_m2,
       vec_m1 = vec_out_m1, vec_m2 = vec_out_m2)
}
reg_results <- get_regression_parameters()

mock_reg_param_estimator <- R6::R6Class("reg_param_estimator", public = list(
  target = target,
  predictors = preds,
  mat_m1 = reg_results$mat_m1,
  mat_m2 = reg_results$mat_m2,
  vec_m1 = reg_results$vec_m1,
  vec_m2 = reg_results$vec_m2,
  parameters = NULL,
  regress = function(modified_2smm, tau, add_var = 0) {
    corr_mat <- self$mat_m1 + self$mat_m2
    corr_mat <- corr_mat + diag(add_var, nrow(corr_mat), ncol(corr_mat))
    corr_vec <- self$vec_m1 + self$vec_m2
    params <- solve(corr_mat, corr_vec)
    names(params) <- colnames(self$predictors)
    self$parameters <- params
  }
))

get_corrected_reg_parms <- function(penalty) {
  wr <- r0 - mean(r0)

  ai11 <- 1
  ai12 <- mean(r0)
  ai13 <- mean(x1)
  ai14 <- mean(r0 * x1)

  ai22 <- mean(r0^2)
  ai23 <- mean(r0 * x1)
  ai24 <- mean(r0^2 * x1)

  ai33 <- mean(x1^2 - sigma[1, 1])
  ai34 <- mean(r0 * x1^2 - r0 * sigma[1, 1])

  ai44 <- mean(r0^2 * x1^2 - r0^2 * sigma[1, 1])

  mat_med <- matrix(c(
    ai11, ai12, ai13, ai14,
    ai12, ai22, ai23, ai24,
    ai13, ai23, ai33, ai34,
    ai14, ai24, ai34, ai44
  ), 4, 4, byrow = TRUE)

  ai1 <- mean(m)
  ai2 <- mean(m * r0)
  ai3 <- mean(m * x1 - sigma[1, 2])
  ai4 <- mean(r0 * m * x1 - r0 * sigma[1, 2])

  vec_med <- c(ai1, ai2, ai3, ai4)

  gamma_med <- solve(mat_med, vec_med)
  wm <- (gamma_med[2] + gamma_med[4] * x1) * wr

  bi11 <- 1
  bi12 <- mean(r0)
  bi13 <- mean(x1)
  bi14 <- mean(m)
  bi15 <- mean(x1 * r0)

  bi21 <- mean(wr)
  bi22 <- mean(wr * r0)
  bi23 <- mean(wr * x1)
  bi24 <- mean(wr * m)
  bi25 <- mean(wr * x1 * r0)

  bi31 <- mean(x1)
  bi32 <- mean(x1 * r0)
  bi33 <- mean(x1^2)
  bi34 <- mean(x1 * m)
  bi35 <- mean(x1^2 * r0)

  bi41 <- mean(wm)
  bi42 <- mean(wm * r0)
  bi43 <- mean((gamma_med[2] * x1 + gamma_med[4] * x1^2) * wr)
  bi44 <- mean((gamma_med[2] * m + gamma_med[4] * x1 * m) * wr)
  bi45 <- mean((gamma_med[2] * x1 + gamma_med[4] * x1^2) * wr * r0)

  bi51 <- mean(x1 * wr)
  bi52 <- mean(x1 * wr * r0)
  bi53 <- mean(x1 * wr * x1)
  bi54 <- mean(x1 * wr * m)
  bi55 <- mean(x1^2 * wr * r0)

  mat_out_m1 <- t(matrix(c(
    bi11, bi21, bi31, bi41, bi51,
    bi12, bi22, bi32, bi42, bi52,
    bi13, bi23, bi33, bi43, bi53,
    bi14, bi24, bi34, bi44, bi54,
    bi15, bi25, bi35, bi45, bi55
  ), 5, 5, byrow = TRUE))
  colnames(mat_out_m1) <- rownames(mat_out_m1) <- pred_names

  ci33 <- mean(- sigma[1, 1])
  ci34 <- mean(- sigma[1, 2])
  ci35 <- mean(- sigma[1, 1] * r0)

  ci43 <- mean(- gamma_med[4] * sigma[1, 1] * wr)
  ci44 <- mean(- gamma_med[4] * sigma[1, 2] * wr)
  ci45 <- mean(-gamma_med[4] * sigma[1, 1] * wr * r0)

  ci55 <- mean(- sigma[1, 1] * wr * r0)

  mat_out_m2 <- t(matrix(c(
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, ci33, ci43, 0,
    0, 0, ci34, ci44, 0,
    0, 0, ci35, ci45, ci55
  ), 5, 5, byrow = TRUE))
  colnames(mat_out_m2) <- rownames(mat_out_m2) <- pred_names

  mat_out <- mat_out_m1 + mat_out_m2

  bi1 <- mean(y)
  bi2 <- mean(wr * y)
  bi3 <- mean(y * x1)
  bi4 <- mean((gamma_med[2] * y + gamma_med[4] * x1 * y) * wr)
  bi5 <- mean(wr * y * x1)

  vec_out_m1 <- c(bi1, bi2, bi3, bi4, bi5)
  names(vec_out_m1) <- pred_names

  ci3 <- mean(- sigma[1, 3])
  ci4 <- mean(- gamma_med[4] * sigma[1, 3] * wr)
  ci5 <- mean(- wr * sigma[1, 3])

  vec_out_m2 <- c(0, 0, ci3, ci4, ci5)
  names(vec_out_m2) <- pred_names

  vec_out <- vec_out_m1 + vec_out_m2

  theta <- setNames(as.numeric(solve(mat_out + penalty * diag(5), vec_out)),
                    pred_names)

  list(theta = theta, mat_m1 = mat_out_m1, mat_m2 = mat_out_m2,
       vec_m1 = vec_out_m1, vec_m2 = vec_out_m2)
}
