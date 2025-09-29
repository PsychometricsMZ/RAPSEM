#' @title R6 Class for Regression Correction based on G-Estimation
#' @description
#' Internal R6 class that implements correction of regression weights and
#' adjustment terms using factor scores, particularly in the context of
#' mediation or structural modeling.
#' @importFrom R6 R6Class
#' @keywords internal
reg_corrector <- R6Class("reg_corrector",
  public = list(
    #' @field corrected_reg_estimator An updated regression estimator object
    #' with applied corrections.
    corrected_reg_estimator = NULL,

    #' @description
    #' Initializes the `reg_corrector` by applying factor score corrections
    #' to the specified regression estimator.
    #' @param reg_est An instance of `reg_param_estimator`, containing
    #'                methods for prediction and regression.
    #' @param reg_model An instance of `reg_model`, containing the
    #'                  mediator predictors.
    #' @param sample A data frame representing the sample data.
    #' @param var_names An instance of `var_names_config`,
    #'                  containing the variable names.
    #' @param factor_scores An instance of `factor_score_estimator`,
    #'                      containing the factor scores and covariances.
    #' @param modified_2smm Boolean. Turns on modified 2SMM estimation.
    #' @param tau Numeric. Modification factor.
    initialize = function(reg_est, reg_model, sample,
                          var_names, factor_scores,
                          modified_2smm, tau) {
      treat_names <- var_names$treatment
      out_pred_names <- colnames(reg_est$predictors)
      if (!any(treat_names %in% out_pred_names)) {
        stop("At least one treatment must be included as a predictor.")
      }
      corr_reg_est <- private$correct_regression(out_pred_names,
                                                 reg_est,
                                                 reg_model,
                                                 sample,
                                                 var_names,
                                                 factor_scores,
                                                 modified_2smm,
                                                 tau)
      self$corrected_reg_estimator <- corr_reg_est
    },

    #' @description
    #' Performs the final regression estimation using the corrected weights and
    #' adjustment terms.
    #' @param modified_2smm Boolean. Turns on modified 2SMM estimation.
    #' @param tau Numeric. Modification factor.
    #' @param add_var Numeric. Small value added to the diagonal of the matrix
    #'                to improve stability (Ridge regularization).
    #' @return Regression coefficients obtained with g-estimation.
    g_estimation = function(modified_2smm, tau, add_var) {
      warnings <- character()
      withCallingHandlers(
        {
          tryCatch(
            {
              result <- self$corrected_reg_estimator$regress(modified_2smm,
                                                             tau,
                                                             add_var)
              return(result)
            },
            error = function(e) {
              stop(paste0("Error during g_estimation: ", e$message),
                   call. = FALSE)
            }
          )
        },
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      if (length(warnings) > 0) {
        warning_msg <- paste0("Warning(s) during g_estimation: ",
                              paste(warnings, collapse = "; "))
        warning(warning_msg, call. = FALSE)
      }
    }
  ),

  private = list(
    #' Correct the regression through g-estimation weights.
    correct_regression = function(out_pred_names, reg_est, reg_model,
                                  sample, var_names, factor_scores,
                                  modified_2smm, tau) {
      out_vars <- cbind(reg_est$predictors, reg_est$target)
      for (pred_name in out_pred_names) {
        if (pred_name %in% names(reg_model$med_preds)) {
          g_calc <- private$compute_mediator_corrections(pred_name,
                                                         reg_model,
                                                         sample,
                                                         out_vars,
                                                         var_names,
                                                         factor_scores,
                                                         modified_2smm,
                                                         tau)
          reg_est$mat_m1[pred_name, ] <- g_calc$vals[-length(g_calc$vals)]
          reg_est$vec_m1[pred_name] <- g_calc$vals[length(g_calc$vals)]
          reg_est$mat_m2[pred_name, ] <- g_calc$corrs[-length(g_calc$corrs)]
          reg_est$vec_m2[pred_name] <- g_calc$corrs[length(g_calc$corrs)]
        } else {
          split_parts <- strsplit(pred_name, ":")
          matched_treatment <- sapply(split_parts, function(parts) {
            matched <- intersect(parts, var_names$treatment)
            if (length(matched) > 0) matched[1] else NA
          })
          if (!is.na(matched_treatment)) {
            g_calc <- private$compute_treatment_corrections(pred_name,
                                                            matched_treatment,
                                                            out_vars,
                                                            factor_scores)
            reg_est$mat_m1[pred_name, ] <- g_calc$vals[-length(g_calc$vals)]
            reg_est$vec_m1[pred_name] <- g_calc$vals[length(g_calc$vals)]
            reg_est$mat_m2[pred_name, ] <- g_calc$corrs[-length(g_calc$corrs)]
            reg_est$vec_m2[pred_name] <- g_calc$corrs[length(g_calc$corrs)]
          }
        }
      }
      return(reg_est)
    },

    #' Applies centering to treatment variables.
    correct_treatment = function(treat_vals) {
      return(treat_vals - mean(treat_vals))
    },

    #' Computes the corrected entries for treatment predictors.
    compute_treatment_corrections = function(pred_name, treat_name,
                                             out_vars, factor_scores) {
      vals <- setNames(rep(0, ncol(out_vars)), colnames(out_vars))
      corrs <- setNames(rep(0, ncol(out_vars)), colnames(out_vars))
      treat_corrected <- private$correct_treatment(out_vars[, treat_name])
      if (pred_name == treat_name) {
        for (pred in colnames(out_vars)) {
          vals[pred] <- mean(treat_corrected * out_vars[, pred])
        }
      } else if (grepl(paste0(treat_name, ":"), pred_name)
                 || grepl(paste0(":", treat_name), pred_name)) {
        other_pred <- gsub(paste0("(:?", treat_name, ":?)"), "", pred_name)
        for (pred in colnames(out_vars)) {
          pred_el <- create_correction_element(
            pred, out_vars, other_pred, out_vars,
            factor_scores$names, factor_scores$sigma
          )
          vals[pred] <- mean(treat_corrected * pred_el$val)
          corrs[pred] <- mean(treat_corrected * pred_el$corr)
        }
      }
      return(list(vals = vals,
                  corrs = corrs))
    },

    #' Computes the corrected entries for mediator predictors.
    compute_mediator_corrections = function(med_name, reg_model, data,
                                            out_vars, name_config,
                                            factor_scores, modified_2smm,
                                            tau) {
      med_pred_names <- reg_model$med_preds[[med_name]]
      med_est <- reg_param_estimator$new(data, name_config, med_name,
                                         factor_scores, med_pred_names)
      gamma <- med_est$regress(modified_2smm, tau)
      vals <- matrix(0, nrow = nrow(out_vars), ncol = ncol(out_vars))
      colnames(vals) <- colnames(out_vars)
      corrs <- matrix(0, nrow = nrow(out_vars), ncol = ncol(out_vars))
      colnames(corrs) <- colnames(out_vars)
      for (i in seq_along(med_pred_names)) {
        med_pred_name <- med_pred_names[i]
        for (treat_name in name_config$treatment) {
          if (treat_name == med_pred_name) {
            treat_corrected <- private$correct_treatment(
              med_est$predictors[, treat_name]
            )
            for (pred in colnames(out_vars)) {
              val <- gamma[i] * treat_corrected * out_vars[, pred]
              vals[, pred] <- vals[, pred] + val
            }
          } else if (grepl(paste0(treat_name, ":"), med_pred_name)
                     || grepl(paste0(":", treat_name), med_pred_name)) {
            treat_corrected <- private$correct_treatment(
              med_est$predictors[, treat_name]
            )
            other_pred <- gsub(paste0("(:?", treat_name, ":?)"), "",
                               med_pred_name)
            for (pred in colnames(out_vars)) {
              pred_el <- create_correction_element(
                pred, out_vars, other_pred, med_est$predictors,
                factor_scores$names, factor_scores$sigma
              )
              val <- gamma[i] * treat_corrected * pred_el$val
              vals[, pred] <- vals[, pred] + val
              corr <- gamma[i] * treat_corrected * pred_el$corr
              corrs[, pred] <- corrs[, pred] + corr
            }
          }
        }
      }
      return(list(vals = apply(vals, 2, mean),
                  corrs = apply(corrs, 2, mean)))
    }
  )
)
