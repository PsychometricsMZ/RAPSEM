#' @title Confirmatory Factor Analysis (CFA) Model Class
#' @description
#' Internal R6 class to extract and organize CFA model indicators.
#' @details
#' Parses lavaan model syntax specific to CFA structure and separates
#' fixed vs. free indicators per factor.
#' @importFrom R6 R6Class
#' @importFrom lavaan lavaanify
#' @keywords internal
cfa_model <- R6Class("cfa_model",
  public = list(
    #' @field lavaan_syntax Lavaan model syntax.
    lavaan_syntax = NULL,
    #' @field fact_names A character vector of factor names
    fact_names = NULL,
    #' @field free_ind A named list of free indicators per factor
    free_ind = NULL,
    #' @field fixed_ind A named list of fixed indicators per factor
    fixed_ind = NULL,

    #' @description
    #' Initialize the cfa_model object.
    #' @param lavaan_model A string in lavaan model syntax.
    initialize = function(lavaan_model) {
      indicators <- private$sort_indicators(lavaan_model)
      self$lavaan_syntax <- lavaan_model
      self$free_ind <- indicators$free
      self$fixed_ind <- indicators$fixed
      self$fact_names <- names(self$free_ind)
    }
  ),

  private = list(
    #' Sorts indicators of each factor in fixed (first) and free
    sort_indicators = function(lavaan_model) {
      parsed <- lavaanify(lavaan_model)
      loading_rows <- parsed[parsed$op == "=~", ]
      latent_vars <- unique(loading_rows$lhs)
      fixed_indicators <- list()
      free_indicators <- list()
      for (factor in latent_vars) {
        rows <- which(parsed$lhs == factor & parsed$op == "=~")
        fixed <- parsed[rows[1], "rhs"]
        fixed_indicators[[factor]] <- fixed
        free <- setdiff(parsed[rows, ]$rhs, fixed)
        free_indicators[[factor]] <- free
      }
      list(
        free = free_indicators,
        fixed = fixed_indicators
      )
    }
  )
)
