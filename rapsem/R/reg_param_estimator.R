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
                                    var_names$intercept)
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
                                            factor_scores$sigma,
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
      params <- private$perform_regression(self$mat_m1,
                                           self$vec_m1,
                                           self$mat_m2,
                                           self$vec_m2,
                                           self$target,
                                           self$var_target_error,
                                           modified_2smm,
                                           tau,
                                           add_var)
      names(params) <- colnames(self$predictors)
      self$parameters <- params
    }
  ),

  private = list(
    #' Perform the regression calculation to estimate parameters.
    #' Solves the system of equations using either the regular
    #' `solve` function or the generalized inverse (`ginv`) from
    #' the `MASS` package if the matrix is singular.
    perform_regression = function(mat_m1, vec_m1, mat_m2, vec_m2, target,
                                  var_err, modified_2smm, tau, add_var) {
      moments <- private$modify_moments(mat_m1, vec_m1, mat_m2, vec_m2, target,
                                        var_err, modified_2smm, tau)
      mat_m <- moments$mat
      mat_m <- mat_m + diag(add_var, nrow(mat_m), ncol(mat_m))
      vec_m <- moments$vec
      params <- tryCatch({
        solve(mat_m, vec_m)
      }, error = function(e) {
        warning("Using ginv due to singular matrix.")
        if (!requireNamespace("MASS", quietly = TRUE)) {
          stop("The 'MASS' package is required but not installed.")
        }
        MASS::ginv(mat_m) %*% vec_m
      })
      return(params)
    },


    modify_moments = function(mat_m1, vec_m1, mat_m2, vec_m2, target,
                              var_err, modified_2smm, tau) {
      factor <- 1
      if (modified_2smm) {
        r1 <- rbind(
          cbind(mat_m1, vec_m1),
          cbind(t(vec_m1), mean(target^2))
        )
        r2 <- rbind(
          cbind(-mat_m2, -vec_m2),
          cbind(t(vec_m2), var_err)
        )
        r1_inv <- MASS::ginv(r1)
        eig <- eigen(r1_inv %*% r2 %*% t(r1_inv), only.values = TRUE)$values
        lambda_max <- max(Re(eig))
        num_obs <- length(target)
        if (lambda_max >= (1 + 1 / num_obs)) {
          factor <- 1 - tau / num_obs
        } else {
          factor <- 1 / lambda_max - 1 / num_obs - tau / num_obs
        }
      }
      mat_m <- mat_m1 + factor * mat_m2
      vec_m <- vec_m1 + factor * vec_m2
      return(list(mat = mat_m, vec = vec_m))
    },

    #' Center values if necessary (if variable is not
    #' a treatment or intercept and not binary).
    center_if_necessary = function(var_name, vals,
                                   treatment_names,
                                   intercept_name) {
      if (!(var_name %in% treatment_names
            || var_name == intercept_name) &&
            (length(unique(vals)) > 2)) {
        return(scale(vals))
      } else {
        return(vals)
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
          values1 <- row_values[, parts[1]]
          values2 <- row_values[, parts[2]]
          predictors[, pred_name] <- values1 * values2
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
