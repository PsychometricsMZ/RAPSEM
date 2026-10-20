#' @title Regression Parameter Estimator Class
#' @description
#' Internal R6 class for performing regression with factor score corrections.
#' The class stores target and predictor variables, calculates factor
#' corrections, and performs the regression to estimate model parameters.
#' @importFrom R6 R6Class
#' @keywords internal
reg_param_estimator <- R6Class("reg_param_estimator",
  public = list(
    #' @field target Numeric vector of the target variable.
    target = NULL,
    #' @field var_target_error Numeric value containing the error variance
    #' the (latent) target variable
    var_target_error = NULL,
    #' @field predictors Matrix of predictor variables.
    predictors = NULL,
    #' @field weights Optional matrix of weights for regression.
    weights = NULL,
    #' @field mat_m1 Regression matrix without corrections.
    mat_m1 = NULL,
    #' @field vec_m1 Regression vector without corrections.
    vec_m1 = NULL,
    #' @field mat_m2 Regression matrix factor score corrections.
    mat_m2 = NULL,
    #' @field vec_m2 Regression vector factor score corrections.
    vec_m2 = NULL,
    #' @field parameters Estimated regression parameters.
    parameters = NULL,
    #' @field sensitivity_param Condition number
    sensitivity_param = NULL,

    #' @description
    #' Initialize the regression parameter estimator.
    #' @param sample Data frame containing the observed data.
    #' @param var_names An instance of `var_names_config`,
    #'                  containing the variable names.
    #' @param target_name String specifying the name of the target variable.
    #' @param factor_scores An instance of `factor_score_estimator`,
    #'                      containing the factor scores and covariances.
    #' @param predictor_names Character vector specifying the names of the
    #'                        predictor variables.
    initialize = function(sample, var_names, target_name,
                          factor_scores, predictor_names) {
      all_values <- cbind(factor_scores$vals, sample)
      all_values[] <- mapply(function(var_name, vals) {
        private$center_if_necessary(var_name, vals,
                                    var_names$treatment,
                                    var_names$intercept,
                                    factor_scores$names,
                                    isTRUE(factor_scores$standardized))
      }, var_name = colnames(all_values), vals = all_values, SIMPLIFY = FALSE)
      target <- matrix(all_values[, target_name], ncol = 1)
      colnames(target) <- target_name
      self$target <- target
      if (target_name %in% factor_scores$names) {
        self$var_target_error <- factor_scores$sigma[target_name, target_name]
      } else {
        self$var_target_error <- 0
      }
      results <- private$prepare_predictors(all_values,
                                            all_values,
                                            target_name,
                                            factor_scores$names,
                                            correction_sigma(factor_scores),
                                            predictor_names)
      self$predictors <- results$preds
      self$mat_m1 <- results$mat1
      self$vec_m1 <- results$vec1
      self$mat_m2 <- results$mat2
      self$vec_m2 <- results$vec2
    },

    #' @description
    #' Perform regression on the target variable using the predictors
    #' and weights
    #' @param modified_2smm Boolean. Turns on modified 2SMM estimation.
    #' @param tau Numeric. Modification factor.
    #' @param add_var Numeric. Small value added to the diagonal of the matrix
    #'                to improve stability (Ridge regularization). Default is 0.
    #' @return Regression coefficients.
    regress = function(modified_2smm, tau, add_var = 0) {
      results <- private$perform_regression(self$mat_m1,
                                            self$vec_m1,
                                            self$mat_m2,
                                            self$vec_m2,
                                            self$target,
                                            self$var_target_error,
                                            modified_2smm,
                                            tau,
                                            add_var)
      self$sensitivity_param <- results$condition_number
      params <- results$params
      names(params) <- colnames(self$predictors)
      self$parameters <- params
      c(params, sensitivity = results$condition_number)
    }
  ),

  private = list(
    #' Perform the regression calculation to estimate parameters.
    #' Solves the system of equations using either the regular
    #' `solve` function or the generalized inverse (`ginv`) from
    #' the `MASS` package if the matrix is singular.
    perform_regression = function(mat_m1, vec_m1, mat_m2, vec_m2, target,
                                  var_err, modified_2smm, tau, add_var) {
      factor <- private$get_factor_to_modify_moments(mat_m1, vec_m1,
                                                     mat_m2, vec_m2,
                                                     target, var_err,
                                                     modified_2smm, tau)
      mat_m <- mat_m1 + factor * mat_m2
      vec_m <- vec_m1 + factor * vec_m2
      mat_m <- mat_m + diag(add_var, nrow(mat_m), ncol(mat_m))
      cond_number <- kappa(mat_m, exact = TRUE)
      params <- MASS::ginv(mat_m) %*% vec_m
      list(
        params = params,
        condition_number = cond_number
      )
    },


    #' Factor that scales the moment correction in the modified 2SMM
    #' estimator of Wall and Amemiya (2000, eqs. 23-24): gamma is the largest
    #' eigenvalue of R1^(-1) R2 (= R1^(-1/2) R2 R1^(-1/2) for symmetric R1; see
    #' largest_relative_eigenvalue for the g-estimation system). The full
    #' correction (1 - tau/n)
    #' is used if 1/gamma >= 1 + 1/n, otherwise it is shrunk to
    #' 1/gamma - 1/n - tau/n. The factor is therefore at most 1 and does not
    #' depend on (non-singular) recodings of the predictors.
    get_factor_to_modify_moments = function(mat_m1, vec_m1,
                                            mat_m2, vec_m2, target,
                                            var_err, modified_2smm, tau) {
      if (!(modified_2smm) || is.null(tau) || is.nan(tau) || tau == 0) {
        1
      } else {
        r1 <- rbind(
          cbind(mat_m1, vec_m1),
          cbind(t(vec_m1), mean(target^2))
        )
        r2 <- rbind(
          cbind(-mat_m2, -vec_m2),
          cbind(-t(vec_m2), var_err)
        )
        gamma_hat <- private$largest_relative_eigenvalue(r1, r2)
        inv_gamma <- if (gamma_hat <= 0) Inf else 1 / gamma_hat
        num_obs <- length(target)
        if (inv_gamma >= (1 + 1 / num_obs)) {
          1 - tau / num_obs
        } else {
          inv_gamma - 1 / num_obs - tau / num_obs
        }
      }
    },

    #' Largest (real part of the) eigenvalue of R1^(-1) R2. For a symmetric
    #' R1 this equals the largest eigenvalue of R1^(-1/2) R2 R1^(-1/2) used by
    #' Wall and Amemiya (2000). In g-estimation the treatment and mediator rows
    #' of R1 are replaced by instrument conditions, so R1 is not symmetric and
    #' only the general form is defined. Falls back to the generalized
    #' inverse if R1 is singular.
    largest_relative_eigenvalue = function(r1, r2) {
      rel <- tryCatch(solve(r1, r2), error = function(e) MASS::ginv(r1) %*% r2)
      max(Re(eigen(rel, only.values = TRUE)$values))
    },

    #' Center values if necessary (if variable is not
    #' a treatment or intercept and not binary). With standardize = TRUE,
    #' observed variables are also scaled to unit variance; factor scores are
    #' only centered, because they are then already in the latent-standardized
    #' metric of their error covariance (see factor_score_estimator).
    center_if_necessary = function(var_name, vals,
                                   treatment_names,
                                   intercept_name,
                                   factor_names = character(),
                                   standardize = FALSE) {
      if (!(var_name %in% treatment_names
            || var_name == intercept_name) &&
            (length(unique(vals)) > 2)) {
        if (standardize && !(var_name %in% factor_names)) {
          scale(vals)
        } else {
          vals - mean(vals)
        }
      } else {
        vals
      }
    },

    #' Prepare the matrix with centered predictors,
    #' including interaction terms.
    prepare_predictors = function(row_values, col_values, target_name,
                                  fact_names, sigma, predictor_names) {
      num_pred <- length(predictor_names)
      num_obs <- nrow(row_values)
      predictors <- matrix(0, nrow = num_obs, ncol = num_pred)
      colnames(predictors) <- predictor_names
      mat_m1 <- matrix(0, nrow = num_pred, ncol = num_pred)
      rownames(mat_m1) <- colnames(mat_m1) <- predictor_names
      mat_m2 <- matrix(0, nrow = num_pred, ncol = num_pred)
      rownames(mat_m2) <- colnames(mat_m2) <- predictor_names
      vec_m1 <- numeric(num_pred)
      names(vec_m1) <- predictor_names
      vec_m2 <- numeric(num_pred)
      names(vec_m2) <- predictor_names
      for (i in seq_len(num_pred)) {
        pred_name <- predictor_names[i]
        if (grepl(":", pred_name)) {
          parts <- strsplit(pred_name, ":")[[1]]
          predictors[, pred_name] <- Reduce(`*`, lapply(parts, function(p) {
            row_values[, p]
          }))
        } else {
          predictors[, pred_name] <- row_values[, pred_name]
        }
        vec_el <- create_correction_element(
          pred_name, row_values, target_name, col_values,
          fact_names, sigma
        )
        vec_m1[i] <- mean(vec_el$val)
        vec_m2[i] <- mean(vec_el$corr)
        for (j in seq_along(predictor_names)) {
          other_pred_name <- predictor_names[j]
          mat_el <- create_correction_element(
            pred_name, row_values, other_pred_name, col_values,
            fact_names, sigma
          )
          mat_m1[i, j] <- mean(mat_el$val)
          mat_m2[i, j] <- mean(mat_el$corr)
        }
      }
      return(list(preds = predictors, mat1 = mat_m1, mat2 = mat_m2,
                  vec1 = vec_m1, vec2 = vec_m2))
    }
  )
)
