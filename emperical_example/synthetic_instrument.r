library(rapsem)
library(dplyr)
library(haven)
library(mice)

set.seed(123)
lambda       <- as.numeric(Sys.getenv("LAMBDA", "-0.8"))  # e.g. LAMBDA=-0.4 Rscript ...
m            <- 20
out_name     <- "post_ct"
med_name     <- "post_ti"
results_dir  <- Sys.getenv("RESULTS_DIR", "emperical_example/results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
out_file     <- file.path(results_dir, sprintf("valid_instrument_results_lambda%.2f.csv", lambda))
dir.create("emperical_example/output", showWarnings = FALSE, recursive = TRUE)

fixed_model <- "
pre_int =~ pre_int1 + pre_int2 + pre_int3 + pre_int4
pre_mot =~ pre_mot1 + pre_mot2 + pre_mot3 + pre_mot4
qalp_tqca =~ qalp_tqca_p1 + qalp_tqca_p2 + qalp_tqca_p3
qalp_tqts =~ qalp_tqts1 + qalp_tqts2 + qalp_tqts3 + qalp_tqts4
post_ct =~ post_ct_p1 + post_ct_p2 + post_ct_p3
post_ti ~ treat + pre_int + pre_mot + qalp_tqca + qalp_tqts + treat:pre_int + treat:pre_mot + treat:qalp_tqts
post_ct ~ treat + post_ti + pre_int + pre_mot + qalp_tqca + qalp_tqts 
"

temp_data <- readRDS("emperical_example/data/tempData.rds")$data

data_full <- read_sav("emperical_example/data/GRAIL1_2024_09_25.sav") %>%
  zap_labels() %>% zap_formats() %>% zap_widths() %>%
  as.data.frame() %>%
  filter(student %in% temp_data$student)

moderator_pattern <- "^pre_mot|^pre_int|^qalp_tqca|^qalp_tqts"

out_items  <- colnames(data_full)[grepl(paste0("^", out_name, "\\d+$"), colnames(data_full))]
out_items  <- out_items[!grepl("_r\\d+$", out_items)]

int_items  <- colnames(data_full)[grepl("^pre_int\\d+$",  colnames(data_full))]
int_items  <- int_items[!grepl("_r\\d+$", int_items)]

mot_items  <- colnames(data_full)[grepl("^pre_mot\\d+$",  colnames(data_full))]
mot_items  <- mot_items[!grepl("_r\\d+$", mot_items)]

tqca_items <- colnames(data_full)[grepl("^qalp_tqca\\d+$", colnames(data_full))]
tqca_items <- tqca_items[!grepl("_r\\d+$", tqca_items)]

tqts_items <- colnames(data_full)[grepl("^qalp_tqts\\d+$", colnames(data_full))]
tqts_items <- tqts_items[!grepl("_r\\d+$", tqts_items)]

all_indicator_items <- c(out_items, int_items, mot_items, tqca_items, tqts_items)

cat("Items collected per construct:\n")
cat(sprintf("  out_items  (%s): %d\n", out_name,    length(out_items)))
cat(sprintf("  int_items  (pre_int):   %d\n",        length(int_items)))
cat(sprintf("  mot_items  (pre_mot):   %d\n",        length(mot_items)))
cat(sprintf("  tqca_items (qalp_tqca): %d\n",        length(tqca_items)))
cat(sprintf("  tqts_items (qalp_tqts): %d\n",        length(tqts_items)))

df_to_impute <- data_full[, all_indicator_items, drop = FALSE]
df_to_impute <- as.data.frame(lapply(df_to_impute,
  function(x) as.numeric(scale(as.numeric(as.character(x))))))

df_to_impute$treat   <- ifelse(as.numeric(as.character(temp_data$colour)) > 0, 1, 0)
df_to_impute$post_ti <- as.numeric(scale(as.numeric(as.character(temp_data[[med_name]]))))

pre_int_scaled <- as.numeric(scale(as.numeric(as.character(temp_data$pre_int))))

df_to_impute$post_ti <- as.numeric(scale(
  df_to_impute$post_ti + lambda * (df_to_impute$treat * pre_int_scaled)
))

pre_mot_scaled   <- as.numeric(scale(rowMeans(df_to_impute[, mot_items,  drop = FALSE], na.rm = TRUE)))
qalp_tqts_scaled <- as.numeric(scale(rowMeans(df_to_impute[, tqts_items, drop = FALSE], na.rm = TRUE)))

df_fstat <- data.frame(
  post_ti   = df_to_impute$post_ti,
  treat     = df_to_impute$treat,
  pre_int   = pre_int_scaled,
  pre_mot   = pre_mot_scaled,
  qalp_tqts = qalp_tqts_scaled
)

lm_restricted <- lm(post_ti ~ treat + pre_int + pre_mot + qalp_tqts, data = df_fstat)
lm_full       <- lm(post_ti ~ treat + pre_int + pre_mot + qalp_tqts +
                      treat:pre_int + treat:pre_mot + treat:qalp_tqts, data = df_fstat)

f_test        <- anova(lm_restricted, lm_full)
interaction_f <- f_test$F[2]
cat(sprintf("Interaction F-statistic (all treat interactions → post_ti): %.4f\n", interaction_f))

make_parcels <- function(df, base, all_items) {
  items <- sample(all_items)
  n     <- length(items)
  g     <- floor(n / 3)
  idx1  <- items[1:g]
  idx2  <- items[(g + 1):(2 * g)]
  idx3  <- items[(2 * g + 1):n]

  df[[paste0(base, "_p1")]] <- rowMeans(df[, idx1, drop = FALSE], na.rm = TRUE)
  df[[paste0(base, "_p2")]] <- rowMeans(df[, idx2, drop = FALSE], na.rm = TRUE)
  df[[paste0(base, "_p3")]] <- rowMeans(df[, idx3, drop = FALSE], na.rm = TRUE)

  df[paste0(base, c("_p1", "_p2", "_p3"))] <- lapply(
    df[paste0(base, c("_p1", "_p2", "_p3"))], function(x) as.numeric(scale(x)))
  df
}

imp_data <- mice(df_to_impute, m = m, method = "pmm", seed = 123)

extract_med_path <- function(fit_i, i, f_stat) {
  sum_med <- fit_i$summary[["gestimation"]][[out_name]]
  if (is.null(sum_med) || !is.data.frame(sum_med)) return(NULL)
  pred_med <- sum_med[rownames(sum_med) == med_name, , drop = FALSE]
  if (nrow(pred_med) == 0) return(NULL)
  se_col    <- if ("Std.Err"  %in% colnames(pred_med)) "Std.Err"  else 2
  lower_pct <- if ("lower_pct" %in% colnames(pred_med)) "lower_pct" else 3
  upper_pct <- if ("upper_pct" %in% colnames(pred_med)) "upper_pct" else 4
  p_col     <- if ("P(>|z|)"  %in% colnames(pred_med)) "P(>|z|)"  else 6
  base_rows <- data.frame(
    imputation    = i,
    path          = rownames(pred_med),
    estimate      = as.numeric(pred_med[, "Estimate"]),
    se            = as.numeric(pred_med[, se_col]),
    lower_pct     = as.numeric(pred_med[, lower_pct]),
    upper_pct     = as.numeric(pred_med[, upper_pct]),
    p_value       = as.numeric(pred_med[, p_col]),
    interaction_f = f_stat,
    lambda        = lambda,
    stringsAsFactors = FALSE
  )
  if ("sensitivity" %in% rownames(sum_med)) {
    sens_obj <- sum_med["sensitivity", , drop = FALSE]
    sens_row <- data.frame(
      imputation    = i,
      path          = "sensitivity",
      estimate      = as.numeric(sens_obj[1, "Estimate"]),
      se            = as.numeric(sens_obj[1, se_col]),
      lower_pct     = as.numeric(sens_obj[, lower_pct]),
      upper_pct     = as.numeric(sens_obj[, upper_pct]),
      p_value       = as.numeric(sens_obj[1, p_col]),
      interaction_f = f_stat,
      lambda        = lambda,
      stringsAsFactors = FALSE
    )
    base_rows <- rbind(base_rows, sens_row)
  }
  base_rows
}

write.csv(
  # same column order as the rows appended by extract_med_path()
  data.frame(imputation = integer(), path = character(),
             estimate = numeric(), se = numeric(),
             lower_pct = numeric(), upper_pct = numeric(), p_value = numeric(),
             interaction_f = numeric(), lambda = numeric()),
  out_file, row.names = FALSE
)

for (i in 1:m) {
  df_complete <- complete(imp_data, i)

  df_complete <- make_parcels(df_complete, out_name,    out_items)
  df_complete <- make_parcels(df_complete, "pre_int",   int_items)
  df_complete <- make_parcels(df_complete, "pre_mot",   mot_items)
  df_complete <- make_parcels(df_complete, "qalp_tqca", tqca_items)
  df_complete <- make_parcels(df_complete, "qalp_tqts", tqts_items)
  fit <- tryCatch(
    rapsem::med_est(fixed_model, df_complete, "treat"),
    error = function(e) { cat(sprintf("Error in imputation %d: %s\n", i, e$message)); NULL }
  )
  if (is.null(fit)) next

  rows <- extract_med_path(fit, i, interaction_f)
  if (!is.null(rows)) {
    write.table(rows, file = out_file, sep = ",",
                append = TRUE, row.names = FALSE,
                col.names = FALSE, quote = FALSE)
  }
  cat(paste0("Completed imputation: ", i, "/", m, "\n"))
}