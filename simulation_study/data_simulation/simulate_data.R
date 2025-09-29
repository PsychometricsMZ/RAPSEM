library(mvtnorm)
library(here)
source(here::here("data_simulation", "name_settings.R"))

simulate_data <- function(coef) {
  # treatment
  treat <- generate_treatment_data(coef)
  # covariates
  cov_scores <- generate_covariate_data(coef)
  # confounder
  conf <- generate_confounder_data(coef)
  # mediator
  med_info <- generate_mediator_data(coef, treat, cov_scores, conf)
  med_scores <- med_info$vals
  # outcome
  out_scores <- generate_outcome_data(coef, treat, med_scores, cov_scores, conf,
                                      med_info$eff)
  # Get observed indicators
  cov_matrix <- do.call(cbind, lapply(seq_along(coef$cov$vals), function(n) {
    get_observed_data(cov_scores[, n, drop = FALSE], coef$cov$vals[[n]])
  }))
  med_indicators <- get_observed_data(med_scores, coef$med)
  out_indicators <- get_observed_data(out_scores, coef$out)
  # Recode treatment
  treat[treat == -1] <- 0
  # Combine all into final data matrix
  data_matrix <- cbind(
    treat,
    cov_matrix,
    med_indicators,
    out_indicators
  )
  return(data_matrix)
}

#' Generates treatment data as independent samples from a Bernoulli
#' distribution with p = 0.5 to simulate random assignment.
generate_treatment_data <- function(coef) {
  if (is.null(coef$num_obs) || coef$num_obs <= 0) {
    stop("Number of observations (num_obs) must be a positive integer.")
  }
  treat <- as.matrix(sample(c(-1, 1), coef$num_obs, replace = TRUE), ncol = 1)
  colnames(treat) <- treatment_name
  return(treat)
}

#' Generates covariate data as samples from a multivariate normal distribution
#' with mean zero and a symmetric correlation matrix where all off-diagonal
#' elements are set to the value of coef$cov$corr.
generate_covariate_data <- function(coef) {
  if (is.null(coef$num_obs) || coef$num_obs <= 0) {
    stop("Number of observations (num_obs) must be a positive integer.")
  }
  num_covs <- length(coef$cov$vals)
  if (length(coef$cov$vals) == 0) {
    stop("Coefficient list must contain at least one covariate.")
  }
  corr <- coef$cov$corr
  if (is.null(corr) || !is.numeric(corr)) {
    stop("Coefficient list must contain a numeric correlation value.")
  }
  covariates <- as.matrix(rmvnorm(coef$num_obs, rep(0, num_covs),
                                  get_cov_corr(coef)))
  colnames(covariates) <- paste0(covariate_name, seq_along(coef$cov$vals))
  return(covariates)
}

get_cov_corr <- function(coef) {
  num_covs <- length(coef$cov$vals)
  corr <- coef$cov$corr
  cov_corr <- diag(num_covs) * (1 - corr) + corr
  return(cov_corr)
}

#' Generates covariate data as samples from a normal distribution
#' with mean zero and standard deviation one.
generate_confounder_data <- function(coef) {
  if (is.null(coef$num_obs) || coef$num_obs <= 0) {
    stop("Number of observations (num_obs) must be a positive integer.")
  }
  conf <- as.matrix(rnorm(coef$num_obs, 0, 1), ncol = 1)
  colnames(conf) <- confounder_name
  return(conf)
}

add_effect <- function(obj, var_name, effect_size, var_vals) {
  if (!is.nan(effect_size) && effect_size != 0) {
    obj$effects[[var_name]] <- effect_size
    obj$data[[var_name]] <- var_vals
    obj$scores <- obj$scores + obj$effects[[var_name]] * var_vals
  }
  return(obj)
}

inter_name <- function(var1_name, var2_name, separator = ":") {
  return(paste0(var1_name, separator, var2_name))
}

add_covariance_entry <- function(cov, val, var1_name, var2_name, var_names) {
  if (var1_name %in% var_names && var2_name  %in% var_names) {
    cov[var1_name, var2_name] <- cov[var2_name, var1_name] <- val
  }
  return(cov)
}

#' Generates mediator data as result of treatment, covariates, and confounder
#' effects, with a noise term added to the final scores.
generate_mediator_data <- function(coef, treat, covs, conf) {
  mediator <- list(
    effects = list(),
    data = list(),
    scores = numeric(coef$num_obs)
  )
  # treatment main effect
  mediator <- add_effect(mediator, treatment_name,
                         coef$treat$eff_on$med$main, treat)
  for (n in seq_len(ncol(covs))) {
    cov <- covs[, n]
    cov_info <- coef$cov$vals[[n]]
    cov_name <- paste0(covariate_name, n)
    # covariate main effect
    mediator <- add_effect(mediator, cov_name,
                           cov_info$eff_on$med$main, cov)
    # covariate-treatment interaction effect
    mediator <- add_effect(mediator, inter_name(cov_name, treatment_name),
                           cov_info$eff_on$med$inter_treat, cov * treat)
  }
  # confounder main effect
  mediator <- add_effect(mediator, confounder_name,
                         coef$conf$eff_on$med$main, conf)
  # confounder-treatment interaction effect
  mediator <- add_effect(mediator, inter_name(confounder_name, treatment_name),
                         coef$conf$eff_on$med$inter_treat, conf * treat)
  # Add noise
  res_var <- get_med_residual_var(coef, mediator$effects)
  med_scores <- mediator$scores + rnorm(coef$num_obs, 0, res_var)
  med_scores <- as.matrix(med_scores, ncol = 1)
  colnames(med_scores) <- mediator_name
  # print(var(med_scores))
  return(list(vals = med_scores, eff = mediator$effects))
}

get_med_residual_var <- function(coef, effects) {
  var_names <- names(effects)
  # Initialize the matrix with all zero entries
  med_cov <- matrix(0, nrow = length(var_names), ncol = length(var_names))
  colnames(med_cov) <- rownames(med_cov) <- var_names
  # Set variances to one
  diag(med_cov) <- 1.0
  # Add covariate correlation
  cov_corr <- coef$cov$corr
  for (i in 1:(length(coef$cov$vals) - 1)) {
    covi_name <- paste0(covariate_name, i)
    for (j in (i + 1):length(coef$cov$vals)) {
      covj_name <- paste0(covariate_name, j)
      # Cov-cov term
      med_cov <- add_covariance_entry(med_cov, cov_corr,
                                      covi_name, covj_name,
                                      var_names)
      # Cov:treat-cov:treat term
      med_cov <- add_covariance_entry(med_cov, cov_corr,
                                      inter_name(covi_name, treatment_name),
                                      inter_name(covj_name, treatment_name),
                                      var_names)
    }
  }
  # print(med_cov)
  # Calculate predicted mediator values
  effects <- as.vector(unlist(effects))
  med_cov_predicted <- t(effects) %*% med_cov %*% effects
  # Return the residual variance
  return(sqrt(1 - med_cov_predicted))
}

#' Generates outcome data as result of treatment, mediator covariates,
#' and confounder effects, with a noise term added to the final scores.
generate_outcome_data <- function(coef, treat, med, covs, conf, med_eff) {
  outcome <- list(
    effects = list(),
    data = list(),
    scores = numeric(coef$num_obs)
  )
  # treatment main effect
  outcome <- add_effect(outcome, treatment_name,
                        coef$treat$eff_on$out$main, treat)
  # mediator main effect
  outcome <- add_effect(outcome, mediator_name,
                        coef$med$eff_on$out$main, med)
  # mediator-treatment interaction effect
  outcome <- add_effect(outcome, inter_name(mediator_name, treatment_name),
                        coef$med$eff_on$out$inter_treat, med * treat)
  for (n in seq_len(ncol(covs))) {
    cov <- covs[, n]
    cov_info <- coef$cov$vals[[n]]
    cov_name <- paste0(covariate_name, n)
    # covariate main effect
    outcome <- add_effect(outcome, cov_name,
                          cov_info$eff_on$out$main, cov)
    # covariate-treatment interaction effect
    outcome <- add_effect(outcome, inter_name(cov_name, treatment_name),
                          cov_info$eff_on$out$inter_treat, cov * treat)
    # covariate-mediator interaction effect
    outcome <- add_effect(outcome, inter_name(cov_name, mediator_name),
                          cov_info$eff_on$out$inter_med, cov * med)
  }
  # confounder main effect
  outcome <- add_effect(outcome, confounder_name,
                        coef$conf$eff_on$out$main, conf)
  # confounder-treatment interaction effect
  outcome <- add_effect(outcome, inter_name(confounder_name, treatment_name),
                        coef$conf$eff_on$out$inter_treat, conf * treat)
  # Add noise (using updated residual variance calculation)
  res_var <- get_out_residual_var(coef, outcome$effects, med_eff)
  out_scores <- outcome$scores + rnorm(coef$num_obs, 0, res_var)
  out_scores <- as.matrix(out_scores, ncol = 1)
  colnames(out_scores) <- outcome_name
  # print(var(out_scores))
  return(out_scores)
}

get_out_residual_var <- function(coef, out_eff, med_eff) {
  var_names <- names(out_eff)
  num_vars <- length(out_eff)
  # Initialize the matrix with all zero entries
  out_cov <- matrix(0, nrow = length(var_names), ncol = length(var_names))
  colnames(out_cov) <- rownames(out_cov) <- var_names
  # Set variances to one
  diag(out_cov) <- 1.0
  # Add covariances of main terms
  cov_corr <- coef$cov$corr
  for (i in 1:(num_vars - 1)) {
    var1_name <- var_names[i]
    for (j in (i + 1):num_vars) {
      var2_name <- var_names[j]
      if (!grepl(":", var1_name) && !grepl(":", var2_name)) {
        if (grepl("cov", var1_name) && grepl("cov", var2_name)) {
          out_cov <- add_covariance_entry(out_cov, cov_corr,
                                          var1_name, var2_name,
                                          var_names)
        } else if (grepl(mediator_name, var1_name)) {
          out_cov <- add_covariance_entry(out_cov, med_eff[[var2_name]],
                                          var1_name, var2_name,
                                          var_names)
        }  else if (grepl(mediator_name, var2_name)) {
          out_cov <- add_covariance_entry(out_cov, med_eff[[var1_name]],
                                          var1_name, var2_name,
                                          var_names)
        }
      }
    }
  }
  # Add covariances of interaction terms
  for (i in 1:(num_vars - 1)) {
    var1_name <- var_names[i]
    if (!grepl(":", var1_name)) next
    parts1 <- unlist(strsplit(var1_name, ":"))
    for (j in (i + 1):num_vars) {
      var2_name <- var_names[j]
      if (!grepl(":", var2_name)) next
      parts2 <- unlist(strsplit(var2_name, ":"))
      cov_el <- (out_cov[parts1[1], parts2[1]]
                 * out_cov[parts1[2], parts2[2]] +
                   out_cov[parts1[1], parts2[2]]
                   * out_cov[parts1[2], parts2[1]])
      out_cov <- add_covariance_entry(out_cov, cov_el,
                                      var1_name, var2_name,
                                      var_names)
    }
  }
  # print(out_cov)
  # Calculate predicted meditor values
  effects <- as.vector(unlist(out_eff))
  out_corr_predicted <- t(effects) %*% out_cov %*% effects
  # Return the residual variance
  return(sqrt(1 - out_corr_predicted))
}

get_observed_data <- function(data_vector, fact_info) {
  if (!fact_info$num_ind > 1) {
    return(data_vector)
  }
  num_obs <- length(data_vector)
  num_ind <- as.numeric(fact_info$num_ind)
  # Randomized item-specific reliabilities
  rel_i <- runif(num_ind, fact_info$rel - 0.1, fact_info$rel + 0.1)
  # Compute adjustment factors for each item (residual variance)
  adjustment_factors <- 1 / rel_i - 1
  # Generate noise: each column gets its own variance
  noise_matrix <- matrix(NA, nrow = num_obs, ncol = num_ind)
  for (i in 1:num_ind) {
    noise_matrix[, i] <- rnorm(num_obs, mean = 0,
                               sd = sqrt(adjustment_factors[i]))
  }
  # Generate random intercepts
  intercepts <- runif(num_ind, min = -1, max = 1)
  # Build observed indicator matrix
  indicator_matrix <- matrix(NA, nrow = num_obs, ncol = num_ind)
  for (i in 1:num_ind) {
    indicator_matrix[, i] <- data_vector + noise_matrix[, i] + intercepts[i]
  }
  name_base <- (if (!is.null(colnames(data_vector))) colnames(data_vector)
                else "var")
  colnames(indicator_matrix) <- paste0(name_base, "_",
                                       indicator_name, 1:num_ind)
  # loadings <- apply(indicator_matrix, 2, function(x) cor(x, data_vector))
  # print(loadings)
  return(indicator_matrix)
}
