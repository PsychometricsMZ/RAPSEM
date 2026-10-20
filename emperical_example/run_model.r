library(rapsem)
library(dplyr)
library(haven)
library(stringr)
library(mice)

set.seed(123)
m <- 20
out_name_summary <- "lo_post_z"
out_name <- "post_ct"
mediator <- "ti"
med_name <- paste0("post_", mediator)
moderator_names <- "^pre_mot|^pre_int|^qalp_tqca|^qalp_tqts"
# RESULTS_DIR (environment variable) keeps runs with different rapsem
# versions apart; default: emperical_example/results
results_dir      <- Sys.getenv("RESULTS_DIR", "emperical_example/results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
out_file         <- file.path(results_dir, "parameters.csv")
out_file_rel     <- file.path(results_dir, "reliability.csv")
dir.create("emperical_example/output", showWarnings = FALSE, recursive = TRUE)

data <- readRDS("emperical_example/data/tempData.rds")
data <- data$data
colnames(data)
df <- data.frame(
  treat = ifelse(as.numeric(as.character(data$colour)) > 0, 1, 0),
  out = scale(data[[out_name_summary]]),
  med = scale(data[[med_name]])
)
moderators <- c(names(data)[grepl(moderator_names, names(data))])
base_moderators <- gsub("_z$", "", moderators)
moderators <- moderators[!duplicated(base_moderators)]
for (mod in moderators) {
  df[[mod]] <- scale(as.numeric(as.character(data[[mod]])))
}

# Initialize results data frames with F-stat columns
results_med <- data.frame(
  moderator = character(),
  coef_int  = numeric(),
  f_stat    = numeric(),
  p_f_test  = numeric(),
  stringsAsFactors = FALSE
)

results_out <- data.frame(
  moderator = character(),
  coef_int  = numeric(),
  f_stat    = numeric(),
  p_f_test  = numeric(),
  stringsAsFactors = FALSE
)

for (mod in moderators) {
  mod_vec <- df[[mod]]
  fit_med_red  <- lm(med ~ treat + mod_vec, data = df)
  fit_med_full <- lm(med ~ treat * mod_vec, data = df)
  f_table_med  <- anova(fit_med_red, fit_med_full)
  results_med <- rbind(results_med, data.frame(
    moderator = mod,
    coef_int  = summary(fit_med_full)$coefficients["treat:mod_vec", "Estimate"],
    f_stat    = f_table_med[2, "F"],
    p_f_test  = f_table_med[2, "Pr(>F)"]
  ))
  fit_out_red  <- lm(out ~ treat + mod_vec + med, data = df)
  fit_out_full <- lm(out ~ treat * mod_vec + med, data = df)
  f_table_out  <- anova(fit_out_red, fit_out_full)
  results_out <- rbind(results_out, data.frame(
    moderator = mod,
    coef_int  = tryCatch(summary(fit_out_full)$coefficients["treat:mod_vec",
                                                            "Estimate"],
                         error = function(e) NA_real_),
    f_stat    = tryCatch(f_table_out[2, "F"],
                         error = function(e) NA_real_),
    p_f_test  = tryCatch(f_table_out[2, "Pr(>F)"], 
                         error = function(e) NA_real_)
  ))
}
results_combined <- merge(
  results_med %>% rename(f_med = f_stat, p_f_med = p_f_test,
                         coef_med = coef_int),
  results_out %>% rename(f_out = f_stat, p_f_out = p_f_test,
                         coef_out = coef_int),
  by = "moderator"
)

sig_mods_med <- results_med %>% filter(p_f_test < 0.05) %>% pull(moderator)
sig_mods_out <- results_out %>% filter(p_f_test < 0.05) %>% pull(moderator)
interacting_confounders <- unique(c(sig_mods_med, sig_mods_out))

if (length(sig_mods_med) > 0) {
  main_formula <- as.formula(paste("med ~ treat +", paste(sig_mods_med, collapse = " + ")))
  int_formula  <- as.formula(paste("med ~ treat * (", paste(sig_mods_med, collapse = " + "), ")"))
  
  fit_red_joint  <- lm(main_formula, data = df)
  fit_full_joint <- lm(int_formula,  data = df)
  joint_med_res  <- anova(fit_red_joint, fit_full_joint)
}

# Initialize an empty character vector to collect lines
latex_lines <- c()

# 1. Header and Preamble
latex_lines <- c(latex_lines, 
  "% --- START LATEX TABLE ---",
  "\\begin{table}[htbp]",
  "  \\centering",
  "  \\caption{Mediator Model: Significant Individual and Joint Moderation Effects}",
  "  \\label{tab:med_moderation}",
  "  \\begin{tabular}{lrrr}",
  "    \\toprule",
  "    \\textbf{Moderator} & \\textbf{Interaction Coef ($B$)} & \\textbf{F-Statistic} & \\textbf{p-value} \\\\",
  "    \\midrule",
  "    \\multicolumn{4}{l}{\\textit{Individual Interaction Tests}} \\\\"
)

# 2. Individual Significant Results
sig_results_table <- results_med %>% filter(moderator %in% sig_mods_med)

if (nrow(sig_results_table) > 0) {
  for (i in 1:nrow(sig_results_table)) {
    mod_name_esc <- gsub("_", "\\\\_", sig_results_table$moderator[i])
    row <- sprintf("    %s & %.3f & %.2f & %.3f \\\\", 
                   mod_name_esc, 
                   sig_results_table$coef_int[i], 
                   sig_results_table$f_stat[i], 
                   sig_results_table$p_f_test[i])
    latex_lines <- c(latex_lines, row)
  }
} else {
  latex_lines <- c(latex_lines, "    \\multicolumn{4}{c}{No individual moderators reached significance ($p < 0.05$)} \\\\")
}


latex_lines <- c(latex_lines,
  "    \\midrule",
  "    \\multicolumn{4}{l}{\\textit{Combined (Omnibus) Interaction Test}} \\\\"
)

if (length(sig_mods_med) > 0) {
  joint_row <- sprintf("    All Significant Combined & --- & %.2f & %.3f \\\\", 
                       joint_med_res[2, "F"], 
                       joint_med_res[2, "Pr(>F)"])
  latex_lines <- c(latex_lines, joint_row)
} else {
  latex_lines <- c(latex_lines, "    All Significant Combined & --- & N/A & N/A \\\\")
}

# 4. Footer and Closing
latex_lines <- c(latex_lines, 
  "    \\bottomrule",
  "    \\addlinespace",
  "    \\multicolumn{4}{l}{\\small \\textit{Note:} $B$ represents the coefficient of the $Treat \\times Moderator$ interaction term.} \\\\",
  "    \\multicolumn{4}{l}{\\small Individual tests use Partial F-tests comparing models with and without the interaction.} \\\\",
  "  \\end{tabular}",
  "\\end{table}",
  "% --- END LATEX TABLE ---"
)

# 5. Final Single Print Call
final_latex_output <- paste(latex_lines, collapse = "\n")
cat(final_latex_output, "\n")


data_full <- read_sav("emperical_example/data/GRAIL1_2024_09_25.sav")
df_new <- data_full %>%
  zap_labels() %>%
  zap_formats() %>%
  zap_widths() %>%
  as.data.frame() %>%
  filter(student %in% data$student) %>%
  mutate(treat = ifelse(as.numeric(as.character(colour)) > 0, 1, 0))

target_prefixes <- paste0(moderator_names, "|^", med_name, "|^", out_name)
all_cols <- colnames(df_new)[grepl(target_prefixes, colnames(df_new))]
all_cols <- all_cols[!grepl("_r\\d+$", all_cols, ignore.case = TRUE)]
base_names <- unique(gsub("\\d+$", "", all_cols))
latent_map <- list()
for (base in base_names) {
  items <- all_cols[grepl(paste0("^", base, "\\d+$"), all_cols)]
  if (length(items) > 5) {
    latent_map[[base]] <- paste0(base, "_p", 1:3)
  } else if (length(items) > 1) {
    latent_map[[base]] <- items
  }
}
measurement_model <- ""
for (lat in names(latent_map)) {
  measurement_model <- paste0(measurement_model, lat, " =~ ",
                              paste(latent_map[[lat]], collapse = " + "), "\n")
}
lat_moderators  <- names(latent_map)[grepl(moderator_names, names(latent_map))]
obs_moderators  <- moderators[!moderators %in% names(latent_map)]
all_moderators  <- c(lat_moderators, obs_moderators)
structural_med_terms <- c("treat", all_moderators)
if (length(sig_mods_med) > 0) {
  structural_med_terms <- c(structural_med_terms,
                            paste0("treat:", sig_mods_med))
}
structural_med_full <- paste(med_name, "~", paste(structural_med_terms,
                                                  collapse = " + "))
structural_out_terms <- c("treat", med_name, all_moderators)
if (length(sig_mods_out) > 0) {
  structural_out_terms <- c(structural_out_terms,
                            paste0("treat:", sig_mods_out))
}
structural_out_full <- paste(out_name, "~", paste(structural_out_terms,
                                                  collapse = " + "))
full_model <- paste(
  measurement_model,
  structural_med_full,
  structural_out_full,
  sep = "\n"
)
cat("=== Generated Model ===\n")
cat(full_model, "\n")
relevant_bases <- c(
  names(latent_map)[names(latent_map) %in% interacting_confounders],
  names(latent_map)[names(latent_map) %in% c(med_name, out_name)]
)
latent_rel <- latent_map[names(latent_map) %in% relevant_bases]
measurement_model <- ""
for (lat in names(latent_rel)) {
  measurement_model <- paste0(measurement_model, lat, " =~ ",
                              paste(latent_rel[[lat]], collapse = " + "), "\n")
}
structural_med_terms <- c("treat", interacting_confounders)
if (length(sig_mods_med) > 0) {
  structural_med_terms <- c(structural_med_terms,
                            paste0("treat:", sig_mods_med))
}
structural_med_reduced <- paste(med_name, "~", paste(structural_med_terms,
                                                     collapse = " + "))
structural_out_terms <- c("treat", med_name, interacting_confounders)
if (length(sig_mods_out) > 0) {
  structural_out_terms <- c(structural_out_terms,
                            paste0("treat:", sig_mods_out))
}
structural_out_reduced <- paste(out_name, "~", paste(structural_out_terms,
                                                     collapse = " + "))
reduced_model <- paste(
  measurement_model,
  structural_med_reduced,
  structural_out_reduced,
  sep = "\n"
)
cat("=== Generated Model ===\n")
cat(reduced_model, "\n")
structural_out_terms <- c("treat", med_name, "pre_int", "qalp_tqts")
structural_out_reduced <- paste(out_name, "~", paste(structural_out_terms,
                                                     collapse = " + "))
reduced_extreme_model <- paste(
  measurement_model,
  structural_med_reduced,
  structural_out_reduced,
  sep = "\n"
)
cat("=== Generated Model ===\n")
cat(reduced_extreme_model, "\n")

existing_items <- all_cols[all_cols %in% colnames(df_new)]
base_vars <- c(med_name, out_name)
existing_bases <- base_vars[base_vars %in% colnames(df_new)]
existing_obs_mods <- obs_moderators[obs_moderators %in% colnames(df_new)]
vars_to_impute <- unique(c(existing_bases, existing_items, existing_obs_mods))
df_to_impute <- df_new[, vars_to_impute]
df_to_impute <- as.data.frame(lapply(df_to_impute,
  function(x) {
    as.numeric(scale(as.numeric(as.character(x))))
  }
))
df_to_impute$treat <- as.numeric(as.character(df_new$treat))
imp_data <- mice(df_to_impute, m = m, method = "pmm", seed = 123)

make_parcels <- function(df, base, all_cols) {
  original_items <- all_cols[grepl(paste0("^", base, "\\d+$"), all_cols)]
  original_items <- original_items[original_items %in% colnames(df)]
  original_items <- sample(original_items)
  n          <- length(original_items)
  group_size <- ceiling(n / 3)
  df[[paste0(base, "_p1")]] <- rowMeans(df[, original_items[1:group_size],
                                           drop = FALSE], na.rm = TRUE)
  df[[paste0(base, "_p2")]] <- rowMeans(
    df[, original_items[(group_size + 1):min(2 * group_size, n)],
       drop = FALSE], na.rm = TRUE
  )
  df[[paste0(base, "_p3")]] <- rowMeans(
    df[, original_items[(2 * group_size + 1):n],
       drop = FALSE], na.rm = TRUE
  )
  p_names        <- paste0(base, "_p", 1:3)
  df[p_names]    <- lapply(df[p_names], function(x) as.numeric(scale(x)))
  df
}

extract_rows <- function(fit_i, mv, meth, out_name, med_name, i) {
  if (!meth %in% names(fit_i$summary)) return(NULL)
  all_rows <- data.frame()
  sum_obj <- fit_i$summary[[meth]][[out_name]]
  if (is.null(sum_obj) || !is.data.frame(sum_obj)) return(NULL)
  predictor_rows <- sum_obj[rownames(sum_obj) != "sensitivity", , drop = FALSE]
  if (nrow(predictor_rows) == 0) return(NULL)
  se_col <- if ("Std.Err" %in% colnames(predictor_rows)) "Std.Err" else 2
  lower_pct <- if ("lower_pct" %in% colnames(predictor_rows)) "lower_pct" else 3
  upper_pct <- if ("upper_pct" %in% colnames(predictor_rows)) "upper_pct" else 4
  p_col  <- if ("P(>|z|)" %in% colnames(predictor_rows)) "P(>|z|)" else 6
  base_rows <- data.frame(
    imputation    = i,
    model_version = mv,
    method        = meth,
    dependent     = out_name,
    path          = rownames(predictor_rows),
    estimate      = as.numeric(predictor_rows[, "Estimate"]),
    se            = as.numeric(predictor_rows[, se_col]),
    lower_pct     = as.numeric(predictor_rows[, lower_pct]),
    upper_pct     = as.numeric(predictor_rows[, upper_pct]),
    p_value       = as.numeric(predictor_rows[, p_col]),
    stringsAsFactors = FALSE
  )
  if (meth == "gestimation") {
    if ("sensitivity" %in% rownames(sum_obj)) {
      sens_obj <- sum_obj["sensitivity", , drop = FALSE]
      sens_row <- data.frame(
        imputation    = i,
        model_version = mv,
        method        = meth,
        dependent     = out_name,
        path          = "sensitivity",
        estimate      = as.numeric(sens_obj[1, "Estimate"]),
        se            = as.numeric(sens_obj[1, se_col]),
        lower_pct     = as.numeric(sens_obj[, lower_pct]),
        upper_pct     = as.numeric(sens_obj[, upper_pct]),
        p_value       = as.numeric(sens_obj[1, p_col]),
        stringsAsFactors = FALSE
      )
      base_rows <- rbind(base_rows, sens_row)
    }
    sum_med <- fit_i$summary[["mediation"]][[med_name]]
    if (!is.null(sum_med) && is.data.frame(sum_med)) {
      pred_med <- sum_med[rownames(sum_med) != "sensitivity", , drop = FALSE]
      if (nrow(pred_med) > 0) {
        se_col_m <- if ("Std.Err" %in% colnames(pred_med)) "Std.Err" else 2
        lower_pct <- if ("lower_pct" %in% colnames(pred_med)) "lower_pct" else 3
        upper_pct <- if ("upper_pct" %in% colnames(pred_med)) "upper_pct" else 4
        p_col_m  <- if ("P(>|z|)" %in% colnames(pred_med)) "P(>|z|)" else 6
        med_rows <- data.frame(
          imputation    = i,
          model_version = mv,
          method        = meth,
          dependent     = med_name,
          path          = rownames(pred_med),
          estimate      = as.numeric(pred_med[, "Estimate"]),
          se            = as.numeric(pred_med[, se_col_m]),
          lower_pct     = as.numeric(pred_med[, lower_pct]),
          upper_pct     = as.numeric(pred_med[, upper_pct]),
          p_value       = as.numeric(pred_med[, p_col_m]),
          stringsAsFactors = FALSE
        )
        base_rows <- rbind(base_rows, med_rows)
      }
    }
  }
  all_rows <- rbind(all_rows, base_rows)
  if (nrow(all_rows) == 0) return(NULL)
  all_rows
}

extract_reliability <- function(fit_i, mv, i) {
  rel <- fit_i$summary$reliability
  if (is.null(rel) || !is.data.frame(rel) || nrow(rel) == 0) return(NULL)
  data.frame(
    imputation    = i,
    model_version = mv,
    factor        = rownames(rel),
    omega         = as.numeric(rel[, "Estimate"]),
    se            = as.numeric(rel[, "Std.Err"]),
    stringsAsFactors = FALSE
  )
}

all_results <- data.frame()
methods_to_save <- c("lavaan", "gestimation", "standard")
model_versions  <- c("full", "reduced", "reduced_extreme")

headers <- data.frame(
  imputation    = integer(),
  model_version = character(),
  method        = character(),
  dependent     = character(),
  path          = character(),
  estimate      = numeric(),
  se            = numeric(),
  lower_pct     = numeric(),
  upper_pct    = numeric(),
  p_value       = numeric()
)
write.csv(headers, out_file, row.names = FALSE)

rel_headers <- data.frame(
  imputation    = integer(),
  model_version = character(),
  factor        = character(),
  omega         = numeric(),
  se            = numeric()
)
write.csv(rel_headers, out_file_rel, row.names = FALSE)

for (i in 1:m) {
  df_complete <- complete(imp_data, i)
  for (base in names(latent_map)) {
    if (any(grepl("_p\\d$", latent_map[[base]]))) {
      df_complete <- make_parcels(df_complete, base, all_cols)
    }
  }
  # blockwise = FALSE: the mediator (post_ti) has only two indicators and is
  # identified only in the joint measurement model.
  fit_rapsem <- function(model) {
    tryCatch(rapsem::med_est(model, df_complete, "treat",
                             group_lavaan = TRUE, blockwise = FALSE),
             error = function(e) {
               message("Imputation ", i, ": ", conditionMessage(e))
               NULL
             })
  }
  fit_full            <- fit_rapsem(full_model)
  fit_reduced         <- fit_rapsem(reduced_model)
  fit_reduced_extreme <- fit_rapsem(reduced_extreme_model)
  if (is.null(fit_full)
      || is.null(fit_reduced)
      || is.null(fit_reduced_extreme)) next
  fit_list      <- list(full = fit_full,
                        reduced = fit_reduced,
                        reduced_extreme = fit_reduced_extreme)
  rows_this_iter <- data.frame()
  rel_this_iter  <- data.frame()
  for (mv in model_versions) {
    for (meth in methods_to_save) {
      new_rows <- extract_rows(fit_list[[mv]], mv, meth, out_name, med_name, i)
      if (!is.null(new_rows)) rows_this_iter <- rbind(rows_this_iter, new_rows)
    }
    sum_med <- fit_list[[mv]]$summary[["mediation"]][[med_name]]
    if (!is.null(sum_med) && is.data.frame(sum_med)) {
      pred_med <- sum_med[rownames(sum_med) != "sensitivity", , drop = FALSE]
      if (nrow(pred_med) > 0) {
        se_col_m <- if ("Std.Err" %in% colnames(pred_med)) "Std.Err" else 2
        lower_pct <- if ("lower_pct" %in% colnames(pred_med)) "lower_pct" else 3
        upper_pct <- if ("upper_pct" %in% colnames(pred_med)) "upper_pct" else 4
        p_col_m  <- if ("P(>|z|)" %in% colnames(pred_med)) "P(>|z|)" else 6
        med_rows <- data.frame(
          imputation    = i,
          model_version = mv,
          method        = "mediation",
          dependent     = med_name,
          path          = rownames(pred_med),
          estimate      = as.numeric(pred_med[, "Estimate"]),
          se            = as.numeric(pred_med[, se_col_m]),
          lower_pct     = as.numeric(pred_med[, lower_pct]),
          upper_pct     = as.numeric(pred_med[, upper_pct]),
          p_value       = as.numeric(pred_med[, p_col_m]),
          stringsAsFactors = FALSE
        )
        rows_this_iter <- rbind(rows_this_iter, med_rows)
      }
    }
    new_rel <- extract_reliability(fit_list[[mv]], mv, i)
    if (!is.null(new_rel)) rel_this_iter <- rbind(rel_this_iter, new_rel)
  }
  if (nrow(rows_this_iter) > 0) {
    write.table(rows_this_iter, file = out_file, sep = ",",
                append = TRUE, row.names = FALSE,
                col.names = FALSE, quote = FALSE)
    all_results <- rbind(all_results, rows_this_iter)
  }
  if (nrow(rel_this_iter) > 0) {
    write.table(rel_this_iter, file = out_file_rel, sep = ",",
                append = TRUE, row.names = FALSE,
                col.names = FALSE, quote = FALSE)
  }
  cat(paste0("Completed Imputation: ", i, "/", m, "\n"))
}