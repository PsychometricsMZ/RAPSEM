library(rapsem)
library(here)
source(here::here("simulation_study/data_simulation",
                  "create_coefficients.R"))
source(here::here("simulation_study/data_simulation",
                  "generate_lavaan_model.R"))
source(here::here("simulation_study/data_simulation",
                  "simulate_data.R"))

num_boot <- 100
args <- commandArgs(trailingOnly = TRUE)

# Expecting arguments in order:
# add_var, num_obs, interaction, true_theta_m, rel, study_num
add_var    <- if (length(args) >= 1) as.numeric(args[1]) else 0.0
num_obs    <- if (length(args) >= 2) as.numeric(args[2]) else 100
interaction   <- if (length(args) >= 3) as.numeric(args[3]) else 0.204
true_theta_m    <- if (length(args) >= 4) as.numeric(args[4]) else 0.0
rel        <- if (length(args) >= 5) as.numeric(args[5]) else 0.75
conf       <- if (length(args) >= 6) as.numeric(args[6]) else 0.0
output_file <- if (length(args) >= 7) args[7] else "results.csv"
conf_med <- if (length(args) >= 8) as.numeric(args[8]) else 0.0
fit_output_file <- if (length(args) >= 9) args[9] else ""

cat("Running with add_var=", add_var,
    " num_obs=", num_obs,
    " interaction=", interaction,
    " true_theta_m=", true_theta_m,
    " rel=", rel,
    " conf=", conf, "\n")

write_to_output_file <- function(out, num_obs, rel, inter,
                                 theta_m, conf,
                                 variable_name = "med",
                                 outcome_name = "out") {
  if (variable_name == "sensitivity") {
    methods <- c("standard", "gestimation")
  } else {
    methods <- c("standard", "gestimation", "lavaan")
  }
  for (method in methods) {
    med_res <- out$summary[[method]][[outcome_name]][variable_name, ]
    if (is.null(med_res) || nrow(med_res) == 0) next
    estimate <- as.numeric(med_res["Estimate"])
    lower_ci <- as.numeric(med_res["2.5%"])
    upper_ci <- as.numeric(med_res["97.5%"])
    se <- as.numeric(med_res["Std.Err"])
    t_val <- estimate / se
    p_val <- 2 * pnorm(-abs(t_val))
    signif <- !(lower_ci <= 0 & upper_ci >= 0)
    row <- data.frame(
      method = method,
      variable = variable_name,
      mean = estimate,
      sd = se,
      lower_pct = lower_ci,
      upper_pct = upper_ci,
      p_value = p_val,
      significance = signif,
      num_obs = num_obs,
      true_theta_m = theta_m,
      interaction = inter,
      confounding = conf,
      conf_med = conf_med,
      reliability = rel,
      stringsAsFactors = FALSE
    )
    write.table(
      row,
      file = output_file,
      sep = ",",
      row.names = FALSE,
      col.names = FALSE,
      append = TRUE
    )
  }
}

write_fit_indices_to_file <- function(out, num_obs, rel, inter, theta_m,
                                      conf) {
  if (fit_output_file == "") return(NULL)
  fit_df <- out$summary$model_fit
  if (is.null(fit_df) || !is.data.frame(fit_df) || nrow(fit_df) == 0) return(NULL)
  fit_df$num_obs      <- num_obs
  fit_df$true_theta_m <- theta_m
  fit_df$interaction  <- inter
  fit_df$confounding  <- conf
  fit_df$conf_med     <- conf_med
  fit_df$reliability  <- rel
  write.table(
    fit_df,
    file      = fit_output_file,
    sep       = ",",
    row.names = FALSE,
    col.names = FALSE,
    append    = TRUE
  )
}

for (replicate in 1:100) {
  seed <- 100000 * replicate
  run <- sprintf("rep=%d n=%d theta_m=%.3f inter=%.3f conf=%.3f conf_med=%.3f rel=%.3f",
                 replicate, num_obs, true_theta_m, interaction, conf, conf_med, rel)
  print(run)
  tryCatch({
    withCallingHandlers({
      coefs <- create_coef(num_obs, rel, interaction, true_theta_m, conf, conf_med)
      lavaan_model <- generate_lavaan_model(coefs)
      sim_data <- simulate_data_with_retry(coefs, seed)
      # named arguments: the 5th positional argument of med_est() is
      # modified_2smm, not add_var. group_lavaan adds the lavaan comparison
      # and its fit indices; the measurement model is fitted blockwise
      # (default of rapsem 0.0.0.2000). Estimates are unstandardized
      # (default; metric of the anchor indicators, like lavaan's est).
      out <- rapsem::med_est(lavaan_model, sim_data, c("treat"),
                             nboot = num_boot, add_var = add_var,
                             group_lavaan = TRUE)
      write_to_output_file(out, num_obs, rel, interaction,
                           true_theta_m, conf)
      write_to_output_file(out, num_obs, rel, interaction,
                           true_theta_m, conf, "sensitivity")
      write_fit_indices_to_file(out, num_obs, rel, interaction, true_theta_m,
                                conf)
    }, warning = function(w) {
      msg <- paste0("WARNING: at ", run, "\n", w$message)
      if (exists("sim_data")) msg <- paste0(msg, "\n", head(sim_data))
      message(msg)
      invokeRestart("muffleWarning")
    })
  }, error = function(e) {
    msg <- paste0("ERROR: at ", run, "\n", e$message)
    if (exists("sim_data")) msg <- paste0(msg, "\n", head(sim_data))
    message(msg)
  })
}
