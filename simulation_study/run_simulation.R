library(rapsem)
library(here)
source(here::here("data_simulation", "create_coefficients.R"))
source(here::here("data_simulation", "generate_lavaan_model.R"))
source(here::here("data_simulation", "simulate_data.R"))

num_boot <- 100
args <- commandArgs(trailingOnly = TRUE)

# Expecting arguments in order:
# add_var, num_obs, beta_mxr, beta_ym, rel, study_num
add_var    <- if (length(args) >= 1) as.numeric(args[1]) else 0.0
num_obs    <- if (length(args) >= 2) as.numeric(args[2]) else 100
beta_mxr   <- if (length(args) >= 3) as.numeric(args[3]) else 0.204
beta_ym    <- if (length(args) >= 4) as.numeric(args[4]) else 0.0
rel        <- if (length(args) >= 5) as.numeric(args[5]) else 0.75
conf       <- if (length(args) >= 6) as.numeric(args[6]) else 0.0
conf_treat <- if (length(args) >= 7) as.numeric(args[7]) else 0.0
output_file <- if (length(args) >= 8) args[8] else "results.csv"

cat("Running with add_var=", add_var,
    " num_obs=", num_obs,
    " beta_mxr=", beta_mxr,
    " beta_ym=", beta_ym,
    " rel=", rel,
    " conf=", conf,
    " conf_treat=", conf_treat, "\n")

write_to_output_file <- function(out, num_obs, rel, beta_mxr,
                                 beta_ym, conf, conf_treat) {
  for (method in c("standard", "gestimation")) {
    med_res <- out$summary[[method]]["med", ]
    estimate <- as.numeric(med_res["Estimate"])
    lower_ci <- as.numeric(med_res["2.5%"])
    upper_ci <- as.numeric(med_res["97.5%"])
    se <- as.numeric(med_res["Std.Err"])
    t_val <- estimate / se
    p_val <- 2 * pnorm(-abs(t_val))
    signif <- !(lower_ci <= 0 & upper_ci >= 0)
    bias <- estimate - beta_ym
    row <- data.frame(
      method = method,
      confounding = conf,
      effect_modification = conf_treat,
      num_obs = num_obs,
      reliability = rel,
      beta_mxr = beta_mxr,
      beta_ym = beta_ym,
      mean = estimate,
      sd = se,
      lower_pct = lower_ci,
      upper_pct = upper_ci,
      p_value = p_val,
      significance = signif,
      bias = bias,
      stringsAsFactors = FALSE
    )
    cat(unlist(row, use.names = FALSE), "\n")
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

for (replicate in 1:100) {
  set.seed(100000 * replicate)
  tryCatch({
    coefs <- create_coef(num_obs, rel, beta_mxr, beta_ym, conf, conf_treat)
    lavaan_model <- generate_lavaan_model(coefs)
    sim_data <- simulate_data(coefs)
    out <- rapsem::med_est(lavaan_model, sim_data, c("treat"),
                           num_boot, add_var)
    write_to_output_file(out, num_obs, rel, beta_mxr, beta_ym, conf, conf_treat)
  }, error = function(e) {
    print(sprintf("Error %s at n=%d rel=%.2f mxr=%.3f ym=%.2f conf=%.2f conf_treat=%.2f rep=%d",
                  e$message, num_obs, rel, beta_mxr, beta_ym, conf, conf_treat,
                  replicate))
  })
}
