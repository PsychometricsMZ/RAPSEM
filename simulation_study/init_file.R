args <- commandArgs(trailingOnly = TRUE)
output_file <- args[1]

column_names <- data.frame(
  method = character(),
  confounding = numeric(),
  effect_modification = numeric(),
  num_obs = integer(),
  reliability = numeric(),
  beta_mxr = numeric(),
  beta_ym = numeric(),
  mean = numeric(),
  sd = numeric(),
  lower_pct = numeric(),
  upper_pct = numeric(),
  p_value = numeric(),
  significance = logical(),
  bias = numeric(),
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
