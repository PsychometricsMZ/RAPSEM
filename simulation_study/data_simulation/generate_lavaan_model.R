library(here)
source(here::here("data_simulation", "name_settings.R"))

generate_lavaan_model <- function(coef) {
  model_lines <- c()

  model_lines <- c(
    model_lines,
    generate_treatment_effects(coef$treat$eff_on),
    generate_covariate_effects(coef$cov$vals),
    generate_mediator_effects(coef$med$eff_on$out),
    generate_measurement_models(coef)
  )

  return(paste(model_lines, collapse = "\n"))
}

generate_treatment_effects <- function(treat) {
  lines <- c()
  if (is_nonzero(treat[[mediator_name]]$main)) {
    lines <- c(lines, paste0(mediator_name, " ~ ", treatment_name))
  }
  if (is_nonzero(treat[[outcome_name]]$main)) {
    lines <- c(lines, paste0(outcome_name, " ~ ", treatment_name))
  }
  return(lines)
}

generate_covariate_effects <- function(covs) {
  lines <- c()
  for (i in seq_along(covs)) {
    cov <- covs[[i]]
    cname <- paste0(covariate_name, i)
    eff <- cov$eff_on

    if (is_nonzero(eff[[mediator_name]]$main)) {
      lines <- c(lines, paste0(mediator_name, " ~ ", cname))
    }
    if (is_nonzero(eff[[outcome_name]]$main)) {
      lines <- c(lines, paste0(outcome_name, " ~ ", cname))
    }
    if (is_nonzero(eff[[mediator_name]]$inter_treat)) {
      lines <- c(lines, paste0(mediator_name, " ~ ",
                               treatment_name, ":", cname))
    }
    if (is_nonzero(eff[[outcome_name]]$inter_treat)) {
      lines <- c(lines, paste0(outcome_name, " ~ ", treatment_name, ":", cname))
    }
    if (is_nonzero(eff[[outcome_name]]$inter_med)) {
      lines <- c(lines, paste0(outcome_name, " ~ ", mediator_name, ":", cname))
    }
  }
  return(lines)
}

generate_mediator_effects <- function(med_eff) {
  lines <- c()
  lines <- c(lines, paste0(outcome_name, " ~ ", mediator_name))
  if (is_nonzero(med_eff$inter_treat)) {
    lines <- c(lines, paste0(outcome_name, " ~ ",
                             treatment_name, ":", mediator_name))
  }
  return(lines)
}

is_nonzero <- function(x) {
  suppressWarnings({
    val <- as.numeric(x)
    return(!is.na(val) && val != 0)
  })
}

generate_measurement_models <- function(coef) {
  lines <- c()
  for (i in seq_along(coef[[covariate_name]]$vals)) {
    cov <- coef[[covariate_name]]$vals[[i]]
    cname <- paste0(covariate_name, i)
    meas <- add_measurement_model(cname, cov$num_ind)
    if (!is.null(meas)) {
      lines <- c(lines, meas)
    }
  }
  med_meas <- add_measurement_model(mediator_name,
                                    coef[[mediator_name]]$num_ind)
  if (!is.null(med_meas)) {
    lines <- c(lines, med_meas)
  }
  out_meas <- add_measurement_model(outcome_name, coef[[outcome_name]]$num_ind)
  if (!is.null(out_meas)) {
    lines <- c(lines, out_meas)
  }
  return(lines)
}

add_measurement_model <- function(varname, num_ind) {
  if (as.numeric(num_ind) > 1) {
    indicators <- paste0(varname, "_", indicator_name, 1:as.numeric(num_ind))
    return(paste0(varname, " =~ ", paste(indicators, collapse = " + ")))
  }
  return(NULL)
}
