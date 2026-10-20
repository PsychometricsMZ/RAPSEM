library(mvtnorm)
library(here)
source(here::here("simulation_study/data_simulation", "name_settings.R"))

simulate_data_with_retry <- function(coef, seed, max_retries = 100) {
  attempt <- 0
  while (attempt < max_retries) {
    current_seed <- seed + attempt
    set.seed(current_seed)
    result <- tryCatch(
      simulate_data(coef),
      error = function(e) {
        if (grepl("variance", conditionMessage(e))) {
          message(sprintf("Seed %d failed: %s Trying seed %d.",
                          current_seed, conditionMessage(e), current_seed + 1))
          return(NULL)
        } else {
          stop(e)
        }
      }
    )
    if (!is.null(result)) return(result)
    attempt <- attempt + 1
    print(paste("attempt:", attempt))
  }
  stop(sprintf("Failed to generate valid data after %d attempts starting from seed %d.",
               max_retries, seed))
}

simulate_data <- function(coef) {
  # treatment
  treat <- generate_treatment_data(coef)
  # covariates
  cov_scores <- generate_covariate_data(coef)
  # confounder
  conf <- generate_confounder_data(coef)
  # mediator
  med_scores <- generate_mediator_data(coef, treat, cov_scores, conf)
  # outcome
  out_scores <- generate_outcome_data(coef, treat, med_scores, cov_scores, conf)
  # Get observed indicators
  cov_matrix <- do.call(cbind, lapply(seq_along(coef$cov$vals), function(n) {
    get_observed_data(cov_scores[, n, drop = FALSE], coef$cov$vals[[n]])
  }))
  med_indicators <- get_observed_data(med_scores, coef$med)
  out_indicators <- get_observed_data(out_scores, coef$out)
  # Combine all into final data matrix
  data_matrix <- cbind(
    treat,
    cov_matrix,
    med_indicators,
    out_indicators
  )
  data_matrix
}

#' Generates treatment data as independent samples from a Bernoulli
#' distribution with p = 0.5 to simulate random assignment.
generate_treatment_data <- function(coef) {
  if (is.null(coef$num_obs) || coef$num_obs <= 0) {
    stop("Number of observations (num_obs) must be a positive integer.")
  }
  treat <- as.matrix(sample(c(-1, 1), coef$num_obs, replace = TRUE), ncol = 1)
  colnames(treat) <- treatment_name
  treat
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
  covariates
}

get_cov_corr <- function(coef) {
  num_covs <- length(coef$cov$vals)
  corr <- coef$cov$corr
  cov_corr <- diag(num_covs) * (1 - corr) + corr
  cov_corr
}

#' Generates covariate data as samples from a normal distribution
#' with mean zero and standard deviation one.
generate_confounder_data <- function(coef) {
  if (is.null(coef$num_obs) || coef$num_obs <= 0) {
    stop("Number of observations (num_obs) must be a positive integer.")
  }
  conf <- as.matrix(rnorm(coef$num_obs, 0, 1), ncol = 1)
  colnames(conf) <- confounder_name
  conf
}

add_effect <- function(obj, var_name, effect_size, var_vals) {
  if (!is.null(effect_size)) {
    if (!is.nan(effect_size) && effect_size != 0) {
      obj$effects[[var_name]] <- effect_size
      obj$data[[var_name]] <- var_vals
      obj$scores <- obj$scores + obj$effects[[var_name]] * var_vals
    }
  }
  obj
}

inter_name <- function(var1_name, var2_name, separator = ":") {
  return(paste0(var1_name, separator, var2_name))
}

add_covariance_entry <- function(cov, val, var1_name, var2_name, var_names) {
  if (var1_name %in% var_names && var2_name  %in% var_names) {
    cov[var1_name, var2_name] <- cov[var2_name, var1_name] <- val
  }
  cov
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
  res_var <- get_residual_var(coef, mediator$effects, "Mediator")
  med_scores <- mediator$scores + rnorm(coef$num_obs, 0, res_var)
  med_scores <- as.matrix(med_scores, ncol = 1)
  colnames(med_scores) <- mediator_name
  # effects of the mediator's causes, needed for its covariances with the
  # other predictors of the outcome (see get_residual_var)
  attr(med_scores, "effects") <- mediator$effects
  #print(var(med_scores))
  vals = med_scores
}

#' Covariance matrix of treatment, covariates, confounder and their
#' treatment interactions (all with mean zero and variance one; treat in
#' {-1, 1}): Cov(cov_i, cov_j) = Cov(cov_i:treat, cov_j:treat) = corr,
#' all other covariances zero.
get_predictor_covariance <- function(coef, var_names) {
  # Initialize the matrix with all zero entries
  covariance <- matrix(0, nrow = length(var_names), ncol = length(var_names))
  colnames(covariance) <- rownames(covariance) <- var_names
  # Set variances to one
  diag(covariance) <- 1.0
  # Add covariate correlation
  cov_corr <- coef$cov$corr
  num_covs <- length(coef$cov$vals)
  for (i in seq_len(num_covs - 1)) {
    covi_name <- paste0(covariate_name, i)
    for (j in (i + 1):num_covs) {
      covj_name <- paste0(covariate_name, j)
      # Cov-cov term
      covariance <- add_covariance_entry(covariance, cov_corr,
                                         covi_name, covj_name,
                                         var_names)
      # Cov:treat-cov:treat term
      covariance <- add_covariance_entry(covariance, cov_corr,
                                         inter_name(covi_name, treatment_name),
                                         inter_name(covj_name, treatment_name),
                                         var_names)
    }
  }
  covariance
}

#' Residual standard deviation such that the latent variable has variance
#' one. In the outcome model, the mediator (variance one) is correlated with
#' its own causes: Cov(med, x) = sum_k b_k Cov(p_k, x) over the mediator's
#' predictors p_k with effects b_k (med_effects). The confounder-mediator
#' product has variance 1 + b_conf^2 and is uncorrelated with the other
#' predictors, given that the confounder enters the mediator only linearly.
get_residual_var <- function(coef, effects, model_name, med_effects = NULL) {
  var_names <- names(effects)
  conf_med <- inter_name(confounder_name, mediator_name)
  unsupported <- setdiff(grep(mediator_name, var_names, value = TRUE),
                         c(mediator_name, conf_med))
  if (length(unsupported) > 0) {
    stop("Residual variance not implemented for: ",
         paste(unsupported, collapse = ", "))
  }
  if (any(c(mediator_name, conf_med) %in% var_names) && is.null(med_effects)) {
    stop("The effects of the mediator's causes are needed.")
  }
  if (!is.null(med_effects[[inter_name(confounder_name, treatment_name)]])) {
    stop("Residual variance not implemented for a confounder-treatment ",
         "interaction in the mediator model.")
  }
  covariance <- get_predictor_covariance(coef, var_names)
  # Mediator-predictor covariances
  if (mediator_name %in% var_names) {
    others <- setdiff(var_names, c(mediator_name, conf_med))
    base <- get_predictor_covariance(coef, union(others, names(med_effects)))
    b_med <- unlist(med_effects)
    cov_med <- as.vector(t(b_med) %*% base[names(b_med), others, drop = FALSE])
    covariance[mediator_name, others] <- cov_med
    covariance[others, mediator_name] <- cov_med
  }
  # Confounder-mediator product
  if (conf_med %in% var_names) {
    b_conf <- if (is.null(med_effects[[confounder_name]])) 0
              else med_effects[[confounder_name]]
    covariance[conf_med, conf_med] <- 1 + b_conf^2
  }
  # Calculate predicted values
  effects <- as.vector(unlist(effects))
  total_variance <- t(effects) %*% covariance %*% effects
  if (as.numeric(total_variance) >= 1) {
    stop(sprintf("%s effects explain %.3f of variance (>= 1).",
                 model_name, total_variance))
  }
  #print(sprintf("%s: effects explain %.3f of variance.",
  #              model_name, total_variance))
  # Return the residual variance
  sqrt(1 - total_variance)
}

#' Generates outcome data as result of treatment, mediator covariates,
#' and confounder effects, with a noise term added to the final scores.
generate_outcome_data <- function(coef, treat, med, covs, conf) {
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
  # confounder-mediator interaction effect
  outcome <- add_effect(outcome, inter_name(confounder_name, mediator_name),
                        coef$conf$eff_on$out$inter_med, conf * med)
  # Add noise (using updated residual variance calculation)
  res_var <- get_residual_var(coef, outcome$effects, "Outcome",
                              attr(med, "effects"))
  out_scores <- outcome$scores + rnorm(coef$num_obs, 0, res_var)
  out_scores <- as.matrix(out_scores, ncol = 1)
  colnames(out_scores) <- outcome_name
  #print(var(out_scores))
  out_scores
}

get_observed_data <- function(data_vector, fact_info) {
  if (!fact_info$num_ind > 1) {
    return(data_vector)
  }
  num_obs <- length(data_vector)
  num_ind <- as.numeric(fact_info$num_ind)
  # Randomized item-specific reliabilities
  rel_i <- pmin(
    pmax(runif(num_ind, fact_info$rel - 0.1, fact_info$rel + 0.1), .01),
    .99
  )
  indicator_matrix <- matrix(NA, nrow = num_obs, ncol = num_ind)
  for (i in seq_len(num_ind)) {
    lambda <- sqrt(rel_i[i])
    error_sd <- sqrt(1 - rel_i[i])
    intercept <- runif(1, -1, 1)
    indicator_matrix[, i] <-
      intercept +
      lambda * data_vector +
      rnorm(num_obs, 0, error_sd)
  }
  if (any(apply(indicator_matrix, 2, var) <= 0)) {
    stop("Indicator variance <= 0")
  }
  name_base <- (if (!is.null(colnames(data_vector))) colnames(data_vector)
                else "var")
  colnames(indicator_matrix) <- paste0(name_base, "_",
                                       indicator_name, 1:num_ind)
  #loadings <- apply(indicator_matrix, 2, function(x) cor(x, data_vector))
  #print(loadings)
  indicator_matrix
}
