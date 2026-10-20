library(dplyr)

params <- read.csv("emperical_example/results/parameters.csv")
m <- max(params$imputation)

sens_data <- params %>%
  filter(method == "gestimation",
         model_version == "full",
         path == "sensitivity")

pooled_robustness <- sens_data %>%
  mutate(
    inv_est = 1 / estimate,
    inv_se  = se / (estimate^2)
  ) %>%
  summarize(
    robustness_mean = mean(inv_est),
    W = mean(inv_se^2),
    B = var(inv_est),
    Total_Var = W + (1 + 1/m) * B,
    pooled_se = sqrt(Total_Var)
  ) %>%
  mutate(
    df      = (m - 1) * (1 + (W / ((1 + 1/m) * B)))^2,
    ci_low  = robustness_mean - qt(0.975, df) * pooled_se,
    ci_high = robustness_mean + qt(0.975, df) * pooled_se
  )

cat(sprintf("Robustness coefficient: %.3f [%.3f, %.3f]\n",
            pooled_robustness$robustness_mean,
            pooled_robustness$ci_low,
            pooled_robustness$ci_high))
