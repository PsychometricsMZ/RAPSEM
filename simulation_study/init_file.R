args <- commandArgs(trailingOnly = TRUE)
output_file <- args[1]

column_names <- data.frame(
  method = character(),
  variable = character(),
  mean = numeric(),
  sd = numeric(),
  lower_pct = numeric(),
  upper_pct = numeric(),
  p_value = numeric(),
  significance = logical(),
  num_obs = integer(),
  true_theta_m = numeric(),
  interaction = numeric(),
  confounding = numeric(),
  conf_med = numeric(),
  reliability = numeric(),
  stringsAsFactors = FALSE
)

cat(names(column_names), "\n")
write.table(
  column_names,
  file = output_file,
  sep = ",",
  row.names = FALSE,
  col.names = TRUE
)
