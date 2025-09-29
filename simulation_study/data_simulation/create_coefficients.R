library(here)
source(here::here("data_simulation", "name_settings.R"))

create_coef <- function(num_obs, rel, beta_mxr, beta_ym, conf, conf_treat) {
  num_indicators <- 3

  result <- list(num_obs = num_obs)

  # treatment
  result[[treatment_name]] <- list(
    eff_on = list(
      med = list(main = 0.3),
      out = list(main = 0.125)
    )
  )

  # covariates
  result[[covariate_name]] <- list(
    corr = 0.2,
    vals = list(
      cov1 = list(
        num_ind = num_indicators,
        rel = rel,
        eff_on = list(
          med = list(main = 0.3, inter_treat = beta_mxr),
          out = list(main = 0.226, inter_treat = 0.0, inter_med = 0.0)
        )
      ),
      cov2 = list(
        num_ind = num_indicators,
        rel = rel,
        eff_on = list(
          med = list(main = 0.3, inter_treat = beta_mxr),
          out = list(main = 0.226, inter_treat = 0.0, inter_med = 0.0)
        )
      )
    )
  )

  # confounder
  result[[confounder_name]] <- list(
    eff_on = list(
      med = list(main = conf, inter_treat = conf_treat),
      out = list(main = conf, inter_treat = conf_treat)
    )
  )

  # mediator
  result[[mediator_name]] <- list(
    num_ind = num_indicators,
    rel = rel,
    eff_on = list(
      out = list(main = beta_ym, inter_treat = 0.0)
    )
  )

  # outcome
  result[[outcome_name]] <- list(
    num_ind = num_indicators,
    rel = rel
  )

  return(result)
}
