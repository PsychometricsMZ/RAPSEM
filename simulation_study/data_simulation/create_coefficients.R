library(here)
source(here::here("simulation_study/data_simulation", "name_settings.R"))

create_coef <- function(num_obs, rel, interaction, true_theta_m, conf,
                        conf_med = 0.0, conf_treat = 0.0) {
  num_indicators <- 3

  result <- list(num_obs = num_obs)

  # treatment
  result[[treatment_name]] <- list(
    eff_on = list(
      out = list(main = 0.224),
      med = list(main = 0.316)
    )
  )

  # covariates
  result[[covariate_name]] <- list(
    corr = 0.3,
    vals = list(
      cov1 = list(
        num_ind = num_indicators,
        rel = rel,
        eff_on = list(
          out = list(main = 0.240, inter_treat = 0, inter_med = 0),
          med = list(main = 0.277, inter_treat = interaction)
        )
      ),
      cov2 = list(
        num_ind = num_indicators,
        rel = rel,
        eff_on = list(
          out = list(main = 0.240, inter_treat = 0, inter_med = 0),
          med = list(main = 0.277, inter_treat = interaction)
        )
      )
    )
  )

  # confounder
  result[[confounder_name]] <- list(
    eff_on = list(
      out = list(main = conf, inter_med = conf_med, inter_treat = conf_treat),
      med = list(main = conf, inter_treat = 0)
    )
  )

  # mediator
  result[[mediator_name]] <- list(
    num_ind = num_indicators,
    rel = rel,
    eff_on = list(
      out = list(main = true_theta_m, inter_treat = 0)
    )
  )

  # outcome
  result[[outcome_name]] <- list(
    num_ind = num_indicators,
    rel = rel
  )

  return(result)
}
