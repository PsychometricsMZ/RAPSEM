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
#' @param tau Numeric. Modification factor. Default is 0.
#' @param add_var Numeric. Small variance value added to improve
#'                stability in g-estimation. Default is 0.
#' @export
med_est <- function(lavaan_txt,
                    data,
                    treat_names,
                    nboot = 100,
                    modified_2smm = TRUE,
                    tau = 5,
                    add_var = 1e-04) {
  data_info <- data_loader$new(data, treat_names)
  models <- model_loader$new(lavaan_txt, data_info$name_config)
  bootstrapped_samples <- data_info$bootstrap(nboot)
  boot_results <- run_bootstrap_reps(
    nboot, bootstrapped_samples, models, data_info, modified_2smm, tau, add_var
  )
  summary_standard <- summarize_results(boot_results$results_standard)
  summary_gestimation <- summarize_results(boot_results$results_gestimation)
  return(list(
    summary = list(
      standard = summary_standard,
      gestimation = summary_gestimation
    )
  ))
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
#' @return A list containing the results
#' @keywords internal
run_bootstrap_reps <- function(nboot, bootstrapped_samples, models, data_info,
                               modified_2smm, tau, add_var) {
  results_standard <- list()
  results_gestimation <- list()

  for (rep in seq_len(nboot)) {
    sample <- bootstrapped_samples[[rep]]
    factor_scores <- factor_score_estimator$new(models$cfa, sample)

    for (out_name in names(models$reg$out_preds)) {
      out_pred_names <- models$reg$out_preds[[out_name]]
      out_est <- reg_param_estimator$new(
        sample,
        data_info$name_config,
        out_name,
        factor_scores,
        out_pred_names
      )
      theta_standard <- out_est$regress(modified_2smm, tau)
      results_standard[[paste0(out_name, "_rep", rep)]] <- theta_standard

      corrected_out_est <- reg_corrector$new(
        out_est,
        models$reg,
        sample,
        data_info$name_config,
        factor_scores,
        modified_2smm,
        tau
      )
      theta_gestimation <- corrected_out_est$g_estimation(modified_2smm, tau,
                                                          add_var)
      results_gestimation[[paste0(out_name, "_rep", rep)]] <- theta_gestimation
    }
  }
  list(
    results_standard = results_standard,
    results_gestimation = results_gestimation
  )
}

#' @title Summarize Bootstrap Estimation Results
#' @description
#' Aggregates bootstrap replicate results by calculating means,
#' standard errors, z-values, p-values, bootstrap percentiles, and
#' adds significance stars for each parameter.
#' @param results_list List of named numeric vectors containing
#' parameter estimates from bootstrap replicates.
#' @return A data.frame with columns Estimate, Std.Err, 2.5% percentile,
#' 97.5% percentile, z value, P(>|z|), and Signif for each parameter.
#' @importFrom stats sd pnorm quantile
#' @keywords internal
summarize_results <- function(results_list) {
  all_params <- unique(unlist(lapply(results_list, names)))
  replicate_matrix <- do.call(rbind, lapply(results_list, function(x) {
    vals <- rep(NA_real_, length(all_params))
    names(vals) <- all_params
    vals[names(x)] <- unlist(x)
    vals
  }))
  est_mean <- colMeans(replicate_matrix, na.rm = TRUE)
  est_se <- apply(replicate_matrix, 2, sd, na.rm = TRUE)
  t_vals <- est_mean / est_se
  p_vals <- 2 * pnorm(-abs(t_vals))
  lower_pct <- apply(replicate_matrix, 2, quantile, probs = 0.025, na.rm = TRUE)
  upper_pct <- apply(replicate_matrix, 2, quantile, probs = 0.975, na.rm = TRUE)
  signif_stars <- function(p) {
    if (is.na(p)) return("")
    else if (p < 0.001) return("***")
    else if (p < 0.01) return("**")
    else if (p < 0.05) return("*")
    else if (p < 0.1) return(".")
    else return("")
  }
  stars <- sapply(p_vals, signif_stars)
  data.frame(
    Estimate = est_mean,
    `Std.Err` = est_se,
    `2.5%` = lower_pct,
    `97.5%` = upper_pct,
    `z value` = t_vals,
    `P(>|z|)` = p_vals,
    Signif = stars,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
