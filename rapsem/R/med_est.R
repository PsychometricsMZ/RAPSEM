#' @title Mediation Estimation for Latent Variables
#' @description
#' Performs mediation estimation using standard factor
#' regression and g-estimation for unbiased estimates in the
#' presence of counfounding with error estimates via bootstrap.
#' @param lavaan_txt lavaan model string
#' @param data input data with column names
#' @param treat_names character vector of treatment names
#' @param nboot number of bootstrap replicates; defaults to 100.
#' @param modified_2smm Boolean. Turns on modified 2SMM estimation.
#'                Default is TRUE.
#' @param tau Numeric. Modification factor of the modified 2SMM estimator
#'            (Wall & Amemiya, 2000). Default is 5.
#' @param add_var Numeric. Small variance value added to improve
#'                stability in g-estimation. Default is 1e-4.
#' @param orthogonal Boolean. Whether to fix the factor covariances to zero
#'                   in the joint measurement model (only used if
#'                   blockwise = FALSE; gives the same estimates as
#'                   blockwise = TRUE). Default is FALSE.
#' @param group_lavaan Boolean. Whether to run group lavaan.
#'                     Default is FALSE.
#' @param blockwise Boolean. Whether to fit the measurement model of each
#'                  latent variable separately (default TRUE). The factor
#'                  scores do not use factor covariances, and separate fits
#'                  keep misspecification in one measurement block from
#'                  affecting the others. Requires at least three indicators
#'                  per latent variable and no residual covariances between
#'                  indicators of different latent variables.
#' @param standardize Boolean. Whether to report standardized coefficients
#'                    (default FALSE). By default, factor scores and their
#'                    error covariance are in the metric of the anchor
#'                    indicators (first indicator of each latent variable,
#'                    loading 1; Wall & Amemiya) and observed variables are
#'                    only centered, like the unstandardized lavaan estimates.
#'                    With TRUE, latent variables are scaled to variance 1
#'                    (estimated latent standard deviation) and observed
#'                    continuous variables to standard deviation 1; the
#'                    treatment and binary variables stay unscaled. The
#'                    lavaan comparison then reports the std.all solution.
#' @export
med_est <- function(lavaan_txt,
                    data,
                    treat_names,
                    nboot = 100,
                    modified_2smm = TRUE,
                    tau = 5,
                    add_var = 1e-04,
                    orthogonal = FALSE,
                    group_lavaan = FALSE,
                    blockwise = TRUE,
                    standardize = FALSE) {
  data_info <- data_loader$new(data, treat_names)
  models <- model_loader$new(lavaan_txt, data_info$name_config)
  warn_unusable_rows(data_info$original_data, models, data_info$name_config)
  bootstrapped_samples <- data_info$bootstrap(nboot)
  boot_results <- run_bootstrap_reps(
    nboot, bootstrapped_samples, models, data_info, modified_2smm,
    tau, add_var, orthogonal, blockwise, standardize
  )
  outcome_predictors <- names(models$reg$out_preds)
  summary <- list(
    standard = summarize_results(boot_results$results_standard,
                                 outcome_predictors),
    gestimation = summarize_results(boot_results$results_gestimation,
                                    outcome_predictors),
    mediation = summarize_results(boot_results$results_mediation,
                                  names(models$reg$med_preds)),
    reliability = summarize_reliability(boot_results$results_reliability)
  )
  if (group_lavaan) {
    fits <- list()
    for (treat_name in treat_names) {
      fit <- lavaan::sem(
        models$group_lavaan,
        data = data,
        group = treat_name,
        group.equal = c("loadings", "intercepts")
      )
      fits[[treat_name]] <- fit
    }
    treated_labels <- vapply(treat_names, function(tr) {
      as.character(max(data[, tr], na.rm = TRUE))
    }, character(1))
    summary[["lavaan"]] <- summarize_lavaan_results(fits, models,
                                                    treated_labels,
                                                    standardize)
    summary[["model_fit"]] <- extract_fit_indices(fits)
  }
  list(summary = summary)
}


#' @title Run Bootstrap Replications for Mediation Estimation
#' @description
#' Performs bootstrap replications of mediation estimation using
#' factor score estimation, standard regression, and g-estimation.
#' @param nboot Integer; number of bootstrap replicates.
#' @param bootstrapped_samples List of bootstrapped data frames.
#' @param models Object containing CFA and regression models.
#' @param data_info Object containing data name configuration.
#' @param modified_2smm Boolean. Turns on modified 2SMM estimation.
#' @param tau Numeric. Modification factor.
#' @param add_var Numeric. Small variance value added to improve
#'                stability in g-estimation.
#' @param orthogonal Boolean. Whether to use orthogonalized factor scores.
#' @param blockwise Boolean. Whether to fit the measurement blocks separately.
#' @param standardize Boolean. Whether to use the standardized metric.
#' @return A list containing the results
#' @keywords internal
run_bootstrap_reps <- function(nboot, bootstrapped_samples, models, data_info,
                               modified_2smm, tau, add_var, orthogonal,
                               blockwise = FALSE, standardize = FALSE) {
  results_standard <- list()
  results_gestimation <- list()
  results_mediation <- list()
  results_reliability <- list()
  failures <- character()
  for (rep in seq_len(nboot)) {
    # A replicate is kept only if all of its estimates succeed.
    rep_results <- tryCatch(
      run_one_bootstrap_rep(rep, bootstrapped_samples[[rep]], models,
                            data_info, modified_2smm, tau, add_var,
                            orthogonal, blockwise, standardize),
      error = function(e) conditionMessage(e)
    )
    if (is.character(rep_results)) {
      failures <- c(failures, rep_results)
      next
    }
    results_reliability[[length(results_reliability) + 1]] <-
      rep_results$reliability
    results_mediation <- c(results_mediation, rep_results$mediation)
    results_standard <- c(results_standard, rep_results$standard)
    results_gestimation <- c(results_gestimation, rep_results$gestimation)
  }
  if (length(failures) == nboot && nboot > 0) {
    stop("All bootstrap replicates failed. First error: ", failures[1])
  }
  if (length(failures) > 0) {
    warning(length(failures), " of ", nboot, " bootstrap replicates failed ",
            "and were skipped. First error: ", failures[1], call. = FALSE)
  }
  list(
    results_standard    = results_standard,
    results_gestimation = results_gestimation,
    results_mediation   = results_mediation,
    results_reliability = results_reliability
  )
}

#' @title Observed Variables of the Structural Model
#' @description
#' Observed (non-latent) variables used in the regression equations,
#' including the parts of interaction terms, excluding the intercept.
#' @param models Object containing CFA and regression models.
#' @param name_config An instance of `var_name_config`.
#' @return Character vector of variable names.
#' @keywords internal
structural_observed_vars <- function(models, name_config) {
  terms <- c(names(models$reg$med_preds), names(models$reg$out_preds),
             unlist(models$reg$med_preds), unlist(models$reg$out_preds))
  parts <- unique(unlist(strsplit(terms, ":", fixed = TRUE)))
  setdiff(parts, c(models$cfa$fact_names, name_config$intercept))
}

#' @title Warn About Unusable Observations
#' @description
#' Counts rows of the original data that cannot be used: a latent variable
#' without any observed indicator, or a missing observed variable of the
#' structural model.
#' @param data Data frame (original data).
#' @param models Object containing CFA and regression models.
#' @param name_config An instance of `var_name_config`.
#' @keywords internal
warn_unusable_rows <- function(data, models, name_config) {
  unusable <- !stats::complete.cases(
    data[, structural_observed_vars(models, name_config), drop = FALSE]
  )
  for (fact in models$cfa$fact_names) {
    inds <- c(models$cfa$fixed_ind[[fact]], models$cfa$free_ind[[fact]])
    unusable <- unusable | rowSums(!is.na(data[, inds, drop = FALSE])) == 0
  }
  if (any(unusable)) {
    warning(sum(unusable), " of ", nrow(data), " observations are not used ",
            "(a latent variable without any observed indicator or a missing ",
            "observed variable of the structural model).", call. = FALSE)
  }
  indicators <- unique(unlist(c(models$cfa$fixed_ind, models$cfa$free_ind)))
  if (anyNA(data[!unusable, indicators])) {
    message("Missing indicator values: the measurement model is fitted with ",
            "full-information ML and factor scores are computed per ",
            "missingness pattern.")
  }
}

#' @title Run One Bootstrap Replicate
#' @description
#' Factor scores, mediator regressions, standard regressions and
#' g-estimation for a single bootstrap sample.
#' @inheritParams run_bootstrap_reps
#' @param rep Integer; replicate number (used in the result names).
#' @param sample Bootstrapped data frame.
#' @return A list with elements reliability, mediation, standard and
#' gestimation.
#' @keywords internal
run_one_bootstrap_rep <- function(rep, sample, models, data_info,
                                  modified_2smm, tau, add_var, orthogonal,
                                  blockwise = FALSE, standardize = FALSE) {
  mediation <- list()
  standard <- list()
  gestimation <- list()
  factor_scores <- factor_score_estimator$new(
    models$cfa, sample, orthogonal,
    keep_vars = structural_observed_vars(models, data_info$name_config),
    blockwise = blockwise,
    standardize = standardize
  )
  # rows without usable factor scores or with missing structural variables
  sample <- sample[factor_scores$rows, , drop = FALSE]
  for (med_name in names(models$reg$med_preds)) {
    med_est <- reg_param_estimator$new(
      sample,
      data_info$name_config,
      med_name,
      factor_scores,
      models$reg$med_preds[[med_name]]
    )
    mediation[[paste0(med_name, "_rep", rep)]] <- med_est$regress(FALSE, tau)
  }
  for (out_name in names(models$reg$out_preds)) {
    out_est <- reg_param_estimator$new(
      sample,
      data_info$name_config,
      out_name,
      factor_scores,
      models$reg$out_preds[[out_name]]
    )
    standard[[paste0(out_name, "_rep", rep)]] <- out_est$regress(FALSE, tau)

    corrected_out_est <- reg_corrector$new(
      out_est,
      models$reg,
      sample,
      data_info$name_config,
      factor_scores,
      modified_2smm,
      tau
    )
    gestimation[[paste0(out_name, "_rep", rep)]] <-
      corrected_out_est$g_estimation(modified_2smm, tau, add_var)
  }
  list(reliability = factor_scores$reliability, mediation = mediation,
       standard = standard, gestimation = gestimation)
}

#' @title Summarize Bootstrap Reliability Estimates
#' @description
#' Aggregates McDonald's omega reliability coefficients across bootstrap
#' replicates. Calculates means, standard errors, z-values, p-values, and
#' bootstrap percentile confidence intervals for each latent factor.
#' @param results_reliability List of named numeric vectors (one per bootstrap
#'   replicate), each containing per-factor McDonald's omega coefficients
#' @return A data.frame with one row per latent factor containing columns
#'   Estimate, Std.Err, 2.5%, 97.5%, z value, and P(>|z|).
#' @importFrom stats sd pnorm quantile
#' @keywords internal
summarize_reliability <- function(results_reliability) {
  replicate_matrix <- do.call(rbind, results_reliability)
  est_mean  <- colMeans(replicate_matrix, na.rm = TRUE)
  est_se    <- apply(replicate_matrix, 2, sd, na.rm = TRUE)
  z_vals    <- est_mean / est_se
  p_vals    <- 2 * pnorm(-abs(z_vals))
  lower_pct <- apply(replicate_matrix, 2, quantile, probs = 0.025, na.rm = TRUE)
  upper_pct <- apply(replicate_matrix, 2, quantile, probs = 0.975, na.rm = TRUE)
  data.frame(
    Estimate  = est_mean,
    `Std.Err` = est_se,
    `2.5%`    = lower_pct,
    `97.5%`   = upper_pct,
    `z value` = z_vals,
    `P(>|z|)` = p_vals,
    check.names    = FALSE,
    stringsAsFactors = FALSE
  )
}

#' @title Summarize Bootstrap Estimation Results by Outcome
#' @description
#' Aggregates bootstrap replicate results separately for each outcome.
#' Calculates means, standard errors, z-values, p-values, and bootstrap
#' percentiles for each parameter within each outcome.
#' @param results_list List of named numeric vectors containing
#' parameter estimates from bootstrap replicates.
#' @param param_names Character vector of parameter names.
#' @return A named list of data.frames (one per outcome) containing
#' Estimate, Std.Err, 2.5%, 97.5%, z value, and P(>|z|).
#' @importFrom stats sd pnorm quantile
#' @keywords internal
summarize_results <- function(results_list, param_names) {
  # results are named "<outcome>_rep<k>"; group them by outcome name
  outcome_of <- sub("_rep[0-9]+$", "", names(results_list))
  split_results <- split(results_list, factor(outcome_of, levels = param_names))
  summarize_one <- function(outcome_results) {
    all_params <- unique(unlist(lapply(outcome_results, names)))
    replicate_matrix <- do.call(rbind, lapply(outcome_results, function(x) {
      vals <- rep(NA_real_, length(all_params))
      names(vals) <- all_params
      vals[names(x)] <- unlist(x)
      vals
    }))
    est_mean <- colMeans(replicate_matrix, na.rm = TRUE)
    est_se <- apply(replicate_matrix, 2, sd, na.rm = TRUE)
    z_vals <- est_mean / est_se
    p_vals <- 2 * pnorm(-abs(z_vals))
    lower_pct <- apply(replicate_matrix, 2, quantile,
                       probs = 0.025, na.rm = TRUE)
    upper_pct <- apply(replicate_matrix, 2, quantile,
                       probs = 0.975, na.rm = TRUE)
    data.frame(
      Estimate = est_mean,
      `Std.Err` = est_se,
      `2.5%` = lower_pct,
      `97.5%` = upper_pct,
      `z value` = z_vals,
      `P(>|z|)` = p_vals,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  lapply(split_results, summarize_one)
}

#' @title Summarize lavaan Estimation Results by Outcome
#' @description
#' Extracts and summarizes parameter estimates from fitted lavaan
#' models separately for each outcome. For each mediator–outcome pair,
#' retrieves estimate, standard error, z-value, p-value, and percentile
#' bootstrap confidence intervals (if available).
#' @param fits Named list of fitted lavaan model objects (one per treatment).
#' @param models Object containing CFA and regression models.
#' @param treated_labels Named character vector giving, per treatment, the
#' value that identifies the treated group. lavaan orders groups by their
#' first appearance in the data, so the group is selected by its label.
#' @param standardize Boolean; if TRUE, the standardized solution (std.all,
#' variances of the treated group) is reported instead of the unstandardized
#' estimates.
#' @return A named list of data.frames (one per outcome). Each data.frame
#' contains columns: Treatment, Mediator, Estimate, Std.Err, 2.5%, 97.5%,
#' z value, and P(>|z|).
#' @importFrom lavaan parameterEstimates lavInspect standardizedSolution
#' @importFrom stats setNames
#' @keywords internal
summarize_lavaan_results <- function(fits, models, treated_labels,
                                     standardize = FALSE) {
  outcome_names <- names(models$reg$out_preds)
  outcome_results <- setNames(vector("list", length(outcome_names)),
                              outcome_names)
  for (outcome_name in outcome_names) {
    outcome_df_list <- list()
    for (treat_name in names(fits)) {
      fit <- fits[[treat_name]]
      group <- match(treated_labels[[treat_name]],
                     lavaan::lavInspect(fit, "group.label"))
      pe <- if (standardize) {
        std <- lavaan::standardizedSolution(fit, type = "std.all",
                                            level = 0.95)
        std$est <- std$est.std
        std
      } else {
        tryCatch(
          lavaan::parameterEstimates(
            fit,
            ci = TRUE,
            level = 0.95,
            boot.ci.type = "perc"
          ),
          error = function(e) {
            lavaan::parameterEstimates(fit)
          }
        )
      }
      target_rows <- pe[
        pe$lhs == outcome_name &
          pe$op == "~" &
          pe$group == group,
      ]
      if (nrow(target_rows) == 0) next
      df <- data.frame(
        Estimate = target_rows$est,
        `Std.Err` = target_rows$se,
        `2.5%` = if ("ci.lower" %in% names(target_rows))
          target_rows$ci.lower else NA_real_,
        `97.5%` = if ("ci.upper" %in% names(target_rows))
          target_rows$ci.upper else NA_real_,
        `z value` = target_rows$z,
        `P(>|z|)` = target_rows$pvalue,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      rownames(df) <- make.unique(target_rows$rhs)
      outcome_df_list[[length(outcome_df_list) + 1]]  <- df
    }
    outcome_results[[outcome_name]] <-
      if (length(outcome_df_list) > 0)
        do.call(rbind, outcome_df_list)
      else
        NULL
  }
  outcome_results
}

#' @title Extract lavaan Model Fit Indices
#' @description
#' Extracts common model fit indices from fitted lavaan model objects.
#' @param fits Named list of fitted lavaan model objects (one per treatment).
#' @return A data.frame containing model fit indices for each treatment model.
#' @importFrom lavaan fitMeasures
#' @keywords internal
extract_fit_indices <- function(fits) {
  fit_list <- lapply(names(fits), function(treat_name) {
    fit <- fits[[treat_name]]
    fm <- lavaan::fitMeasures(
      fit,
      c(
        "chisq",
        "df",
        "pvalue",
        "cfi",
        "tli",
        "rmsea",
        "rmsea.ci.lower",
        "rmsea.ci.upper",
        "srmr",
        "aic",
        "bic"
      )
    )
    data.frame(
      Treatment = treat_name,
      ChiSq = fm["chisq"],
      DF = fm["df"],
      `P(ChiSq)` = fm["pvalue"],
      CFI = fm["cfi"],
      TLI = fm["tli"],
      RMSEA = fm["rmsea"],
      `RMSEA 2.5%` = fm["rmsea.ci.lower"],
      `RMSEA 97.5%` = fm["rmsea.ci.upper"],
      SRMR = fm["srmr"],
      AIC = fm["aic"],
      BIC = fm["bic"],
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, fit_list)
}
