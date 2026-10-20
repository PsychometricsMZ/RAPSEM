library(dplyr)
library(tidyr)

rel_file <- "emperical_example/results/reliability.csv"
out_tex  <- "emperical_example/tables/reliability.tex"

pool_rubin <- function(estimates, ses) {
  m         <- length(estimates)
  theta     <- mean(estimates, na.rm = TRUE)
  within_v  <- mean(ses^2,     na.rm = TRUE)
  between_v <- var(estimates,  na.rm = TRUE)
  total_v   <- within_v + (1 + 1 / m) * between_v
  se        <- sqrt(total_v)
  z         <- theta / se
  p         <- 2 * pnorm(-abs(z))
  data.frame(estimate = theta, se = se, z = z, p_value = p)
}

rel_raw <- read.csv(rel_file, stringsAsFactors = FALSE)

rel_pooled <- rel_raw %>%
  filter(model_version == "full") %>%
  group_by(factor) %>%
  summarise(pool_rubin(omega, se), .groups = "drop") %>%
  rename(omega = estimate)

build_latex_table <- function(rel_pooled) {
  lines <- c(
    "\\begin{table}[htbp]",
    "  \\centering",
    "  \\caption{Pooled McDonald's $\\omega$ Reliability Coefficients (Full Model)}",
    "  \\label{tab:reliability}",
    "  \\begin{tabular}{lrr}",
    "    \\toprule",
    "    Factor & $\\hat{\\omega}$ & SE \\\\",
    "    \\midrule"
  )

  for (i in seq_len(nrow(rel_pooled))) {
    row   <- rel_pooled[i, ]
    label <- gsub("_", "\\\\_", row$factor)
    lines <- c(lines, sprintf(" %s & %.2f & %.2f \\\\",
                              label, row$omega, row$se))
  }

  lines <- c(lines,
    "    \\bottomrule",
    "  \\end{tabular}",
    "\\end{table}"
  )
  paste(lines, collapse = "\n")
}

tex <- build_latex_table(rel_pooled)
writeLines(tex, out_tex)
