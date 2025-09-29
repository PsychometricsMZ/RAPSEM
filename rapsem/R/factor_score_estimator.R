#' @title Factor Score Estimator Class
#' @description
#' Internal R6 class for computing factor scores from lavaan cfa output
#' according to the method of moments technique of Wall and Amemiya (2003).
#' @details
#' Takes a `cfa_model` object and a data sample as input, uses the lavaan
#' `cfa` function with MLR estimation to fit the CFA model to the data
#' and computes the factor scores as well as the the error covariance matrix
#' based on extracted intercepts, loadings, and residual variances with
#' a method of  moments approach.
#' @importFrom R6 R6Class
#' @importFrom lavaan cfa inspect coef
#' @keywords internal
factor_score_estimator <- R6Class("factor_score_estimator",
  public = list(
    #' @field names Character vector of latent factor names.
    names = NULL,
    #' @field vals Matrix of computed factor scores for each observation.
    vals = NULL,
    #' @field sigma Covariance matrix of factor score estimation errors.
    sigma = NULL,

    #' @description
    #' Initializes a `factor_score_estimator` object and computes factor scores
    #' from a CFA model and data sample.
    #' @param cfa_model An instance of `cfa_model` defining the
    #'                  factor structure.
    #' @param data_sample A data frame containing observed indicator values.
    initialize = function(cfa_model, data_sample) {
      tryCatch({
        self$names <- cfa_model$fact_names
        results <- private$calculate_factor_scores(cfa_model, data_sample)
        self$vals <- results$scores
        self$sigma <- results$sigma_ee
      }, error = function(e) {
        stop(paste("Error initializing factor score estimator: ",
                   conditionMessage(e)))
      })
    }
  ),

  private = list(
    #' Core logic for estimating factor scores from a CFA model and data sample.
    calculate_factor_scores = function(cfa_model, data) {
      free_indicators <- unlist(cfa_model$free_ind, use.names = FALSE)
      fixed_indicators <- unlist(cfa_model$fixed_ind, use.names = FALSE)
      all_indicators <- c(free_indicators, fixed_indicators)
      z_values <- data[, all_indicators]
      cfa_coef <- private$get_cfa_coef(cfa_model, z_values)
      beta0hat <- matrix(cfa_coef[paste0(free_indicators, "~1")],
                         nrow = length(free_indicators), ncol = 1)
      beta1hat <- private$build_loading_matrix(cfa_coef, cfa_model,
                                               free_indicators)
      psi_hat <- private$build_residual_covariance_matrix(cfa_coef,
                                                          all_indicators)
      num_facts <- length(cfa_model$fact_names)
      b1_mod1 <- cbind(matrix(0, num_facts, length(free_indicators)),
                       diag(num_facts))
      b1_mod2 <- rbind(diag(length(free_indicators)), -1 * t(beta1hat))
      gamma <- (b1_mod1 %*% psi_hat %*% b1_mod2 %*%
                  solve(t(b1_mod2) %*% psi_hat %*% b1_mod2))
      sigma_ee <- (cbind(-gamma, diag(num_facts) + gamma %*% beta1hat)
                   %*% psi_hat %*% t(b1_mod1))
      num_obs <- nrow(data)
      fact_scores <- matrix(0, nrow = num_obs, ncol = num_facts)
      for (i in 1:num_obs) {
        centered <- t(as.matrix(z_values[i, ])) - rbind(beta0hat,
                                                        matrix(0, num_facts, 1))
        fact_scores[i, ] <- (cbind(-gamma, diag(num_facts) + gamma %*% beta1hat)
                             %*% centered)
      }
      colnames(fact_scores) <- cfa_model$fact_names
      colnames(sigma_ee) <- rownames(sigma_ee) <- cfa_model$fact_names
      return(list(
        scores = fact_scores,
        sigma_ee = sigma_ee
      ))
    },

    #' Fit cfa model, check if it converged and return coefficients
    get_cfa_coef = function(model, indicator_matrix) {
      cfa <- withCallingHandlers(
        cfa(model$lavaan_syntax, data = indicator_matrix,
            meanstructure = TRUE, estimator = "MLR",
            missing = "ML", se = "none", std.lv = TRUE),
        warning = function(w) {
          if (grepl("Model estimation FAILED", conditionMessage(w))) {
            invokeRestart("muffleWarning")
          }
        }
      )
      if (!inspect(cfa, "converged")) {
        stop("CFA model did not converge; cannot compute factor scores.")
      } else {
        return(coef(cfa))
      }
    },

    #' Extract factor loading matrix from lavaan fit
    build_loading_matrix = function(cfa_coef, cfa_model, free_indicators) {
      factor_loadings <- matrix(0, nrow = length(free_indicators),
                                ncol = length(cfa_model$fact_names))
      row_start <- 1
      for (j in seq_along(cfa_model$fact_names)) {
        fact <- cfa_model$fact_names[j]
        indicators <- cfa_model$free_ind[[fact]]
        loading_names <- paste0(fact, "=~", indicators)
        num_indicators <- length(indicators)
        row_end <- row_start + num_indicators - 1
        factor_loadings[row_start:row_end, j] <- cfa_coef[loading_names]
        row_start <- row_end + 1
      }
      return(factor_loadings)
    },

    #' Extract residual covariance matrix from lavaan fit
    build_residual_covariance_matrix = function(cfa_coef, indicators) {
      cov_labels <- outer(indicators, indicators,
        FUN = Vectorize(function(i, j) {
          if (i <= j) paste0(i, "~~", j) else paste0(j, "~~", i)
        })
      )
      cov_values <- matrix(0, nrow = length(indicators),
                           ncol = length(indicators))
      for (i in seq_along(indicators)) {
        for (j in seq_along(indicators)) {
          label <- cov_labels[i, j]
          if (label %in% names(cfa_coef)) {
            cov_values[i, j] <- cfa_coef[label]
          }
        }
      }
      cov_values <- (cov_values + t(cov_values)) / 2
      colnames(cov_values) <- rownames(cov_values) <- indicators
      return(cov_values)
    }
  )
)
