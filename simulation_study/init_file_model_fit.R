args <- commandArgs(trailingOnly = TRUE)
fit_output_file <- args[1]

fit_column_names <- data.frame(
  Treat_Var = character(),
  ChiSq = numeric(),
  DF = numeric(),
  `P(ChiSq)` = numeric(),
  CFI = numeric(),
  TLI = numeric(),
  RMSEA = numeric(),
  `RMSEA 2.5%` = numeric(),
  `RMSEA 97.5%` = numeric(),
  SRMR = numeric(),
  AIC = numeric(),
  BIC = numeric(),
  num_obs = integer(),
  true_theta_m = numeric(),
  interaction = numeric(),
  confounding = numeric(),
  conf_med = numeric(),
  reliability = numeric(),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(names(fit_column_names), "\n")

write.table(
  fit_column_names,
  file = fit_output_file,
  sep = ",",
  row.names = FALSE,
  col.names = TRUE
)