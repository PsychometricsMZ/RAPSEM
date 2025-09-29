#' @title Factor Score Correction Creation
#' @description
#' Internal function for calculating corrected elements in latent
#' regression models based on main effects and interaction terms.
#' @param var_1 The name of the first covariate
#' (either main effect or part of an interaction).
#' @param vals_1 The predictor matrix of the first term.
#' @param var_2 The name of the second covariate
#' (either main effect or part of an interaction).
#' @param vals_2 The predictor matrix of the second term.
#' @param factor_names A character vector of factor (latent) variable names.
#' @param sigma_ee A matrix of factor score covariances,
#' indexed by the factor variable names.
#' @return List containing the value and the correction of the element.
#' @keywords internal
create_correction_element <- function(var_1, vals_1, var_2, vals_2,
                                      factor_names, sigma_ee) {
  if (!grepl(":", var_1) && !grepl(":", var_2)) {
    return(handle_main_effect(var_1, vals_1, var_2, vals_2,
                              factor_names, sigma_ee))
  }
  if (grepl(":", var_1) && !grepl(":", var_2)) {
    return(handle_interaction(var_1, vals_1, var_2, vals_2,
                              factor_names, sigma_ee))
  }
  if (!grepl(":", var_1) && grepl(":", var_2)) {
    return(handle_interaction(var_2, vals_2, var_1, vals_1,
                              factor_names, sigma_ee))
  }
  return(handle_double_interaction(var_1, vals_1, var_2, vals_2,
                                   factor_names, sigma_ee))
}

#' @title Handle Main Effect
#' @description
#' Calculates the corrected vector for a pair of main effect variables.
#' Returns the mean of their product, which is corrected with the corresponding
#' value from the sigma_ee matrix if both variables are found in factor_names.
#' @param var_1 The name of the first main effect variable.
#' @param vals_1 The predictor matrix of the first term.
#' @param var_2 The name of the second main effect variable.
#' @param vals_2 The predictor matrix of the second term.
#' @param factor_names A character vector of factor (latent) variable names.
#' @param sigma_ee A matrix of factor score covariances.
#' @return List containing the value and the correction of the
#' main effects.
#' @keywords internal
handle_main_effect <- function(var_1, vals_1, var_2, vals_2,
                               factor_names, sigma_ee) {
  if (var_1 %in% factor_names && var_2 %in% factor_names) {
    return(
      list(val = vals_1[, var_1] * vals_2[, var_2],
           corr = - sigma_ee[var_1, var_2])
    )
  } else {
    return(
      list(val = vals_1[, var_1] * vals_2[, var_2],
           corr = 0)
    )
  }
}

#' @title Extract Interaction Parts
#' @description
#' Extracts the latent and observed components from a given interaction term.
#' An interaction term can only contain one latent variable.
#' @param inter The interaction term to extract parts from.
#' @param factor_names A character vector of factor (latent) variable names.
#' @return A list with two components: `latent` (latent variables) and
#' `observed` (observed variables).
#' @keywords internal
extract_interaction_parts <- function(inter, factor_names) {
  parts <- unlist(strsplit(inter, ":"))
  latent <- parts[parts %in% factor_names]
  observed <- setdiff(parts, latent)
  if (length(latent) > 1 || length(observed) == 0) {
    stop(paste("Interaction'", inter,
               "can only contain one latent variable."))
  }
  return(list(latent = latent, observed = observed))
}

#' @title Handle Interaction
#' @description
#' Calculates the corrected vector for a main effect and an
#' interaction term. Returns the mean of their product, which is corrected
#' with the corresponding value from the sigma_ee matrix if both variables
#' contain a name found in factor_names.
#' @param inter The interaction term to handle.
#' @param vals_inter The predictor matrix of the interaction term.
#' @param single The single variable (either latent or observed).
#' @param vals_single The predictor matrix of the single term.
#' @param factor_names A character vector of factor (latent) variable names.
#' @param sigma_ee A matrix of factor score covariances.
#' @return List containing the value and the correction of the
#' interaction.
#' @keywords internal
handle_interaction <- function(inter, vals_inter, single, vals_single,
                               factor_names, sigma_ee) {
  parts <- extract_interaction_parts(inter, factor_names)
  product_observed <- apply(vals_inter[, parts$observed, drop = FALSE], 1, prod)
  if (single %in% factor_names && length(parts$latent) == 1) {
    return(
      list(val = product_observed * vals_inter[, parts$latent]
           * vals_single[, single],
           corr = - product_observed * sigma_ee[parts$latent, single])
    )
  } else if (length(parts$latent) == 1) {
    return(
      list(val = product_observed * vals_inter[, parts$latent]
           * vals_single[, single],
           corr = 0)
    )
  } else {
    return(
      list(val = product_observed * vals_single[, single],
           corr = 0)
    )
  }
}

#' @title Handle Double Interaction
#' @description
#' Calculates the corrected vector for a pair of interaction terms.
#' Returns the mean of their product, which is corrected with the
#' corresponding value from the sigma_ee matrix if both variables
#' contain a name found in factor_names.
#' @param var_1 The name of the first interaction term.
#' @param vals_1 The predictor matrix of the first term.
#' @param var_2 The name of the second interaction term.
#' @param vals_2 The predictor matrix of the second term.
#' @param factor_names A character vector of factor (latent) variable names.
#' @param sigma_ee A matrix of factor score covariances.
#' @return List containing the value and the correction of the
#' double interaction.
#' @keywords internal
handle_double_interaction <- function(var_1, vals_1,
                                      var_2, vals_2,
                                      factor_names, sigma_ee) {
  parts_1 <- extract_interaction_parts(var_1, factor_names)
  product_observed_1 <- apply(vals_1[, parts_1$observed, drop = FALSE], 1, prod)
  parts_2 <- extract_interaction_parts(var_2, factor_names)
  product_observed_2 <- apply(vals_2[, parts_2$observed, drop = FALSE], 1, prod)
  if (length(parts_1$latent) == 1 && length(parts_2$latent) == 1) {
    return(
      list(val = product_observed_1 * product_observed_2
           * vals_1[, parts_1$latent] * vals_2[, parts_2$latent],
           corr = - product_observed_1 * product_observed_2
           * sigma_ee[parts_1$latent, parts_2$latent])
    )
  } else if (length(parts_1$latent) == 1) {
    return(
      list(val = product_observed_1 * vals_1[, parts_1$latent]
           * product_observed_2,
           corr = 0)
    )
  } else if (length(parts_2$latent) == 1) {
    return(
      list(val = product_observed_1 * product_observed_2
           * vals_2[, parts_2$latent],
           corr = 0)
    )
  } else {
    return(
      list(val = product_observed_1 * product_observed_2,
           corr = 0)
    )
  }
}
