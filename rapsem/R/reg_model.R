#' @title Regression Model Class
#' @description
#' Internal R6 class to parse and organize regression model predictors.
#' @details
#' Processes lavaan regression statements to identify mediators and outcomes.
#' Adds intercepts where not explicitly omitted.
#' @importFrom R6 R6Class
#' @keywords internal
reg_model <- R6Class("reg_model",
  public = list(
    #' @field text A character vector of regression model lines.
    text = NULL,
    #' @field med_preds A named list of mediator variables and their predictors.
    med_preds = NULL,
    #' @field out_preds A named list of outcome variables and their predictors.
    out_preds = NULL,

    #' @description
    #' Initialize the regression model object.
    #' @param reg_model_lines A character vector of regression syntax lines.
    #' @param intercept_name Name used to represent the intercept term.
    initialize = function(reg_model_lines, intercept_name) {
      self$text <- reg_model_lines
      predictors <- private$get_predictors(reg_model_lines, intercept_name)
      self$med_preds <- predictors$mediators
      self$out_preds <- predictors$outcomes
    }
  ),

  private = list(
    #' Extracts predictors for each regression equation and identifies
    #' mediators and outcomes based on appearance in left- and right-hand sides.
    get_predictors = function(reg_text_lines, intercept_name) {
      lhs_vars <- c()
      rhs_vars <- c()
      regression_info <- list()
      known_covariates <- c()
      for (line in reg_text_lines) {
        if (private$is_regression_line(line)) {
          parsed <- private$parse_regression_line(line, intercept_name,
                                                  known_covariates)
          outcome <- parsed$outcome
          covariates <- parsed$covariates

          lhs_vars <- c(lhs_vars, outcome)
          rhs_vars <- c(rhs_vars, covariates)
          known_covariates <- unique(c(known_covariates, covariates))

          regression_info[[outcome]] <- private$update_regression_info(
            regression_info[[outcome]], covariates
          )
        }
      }
      output <- private$classify_variables(lhs_vars, rhs_vars, regression_info)
      return(output)
    },

    #' Determines whether a regression line specifies a zero intercept.
    has_zero_intercept = function(line) {
      grepl("^\\s*\\w+\\s*~\\s*0\\s*(\\+\\s*\\w+.*|$)", line)
    },

    #' Checks if a line is a valid regression line.
    is_regression_line = function(line) {
      grepl("~", line) && !grepl("=~", line)
    },

    #' Parses a regression line into its outcome and covariates,
    #' handling the intercept and interaction terms.
    parse_regression_line = function(line, intercept_name, known_covariates) {
      parts <- strsplit(line, "~")[[1]]
      outcome <- trimws(parts[1])
      covariates <- unlist(strsplit(parts[2], "\\+"))
      covariates <- trimws(covariates)

      if (private$has_zero_intercept(line)) {
        covariates <- covariates[covariates != "0"]
      } else {
        covariates <- c(intercept_name, covariates)
      }

      expanded_covariates <- covariates
      for (term in covariates) {
        if (grepl(":", term)) {
          parts <- unlist(strsplit(term, ":"))
          for (p in parts) {
            if (!(p %in% covariates) && !(p %in% known_covariates)) {
              expanded_covariates <- c(expanded_covariates, p)
              warning(paste("Main effect", p, "added as predictor",
                            "because interaction", term, "was detected",
                            "in outcome", outcome))
            }
          }
        }
      }
      expanded_covariates <- unique(expanded_covariates)
      list(outcome = outcome, covariates = expanded_covariates)
    },

    #' Updates the regression info list by merging new covariates.
    update_regression_info = function(existing, covariates) {
      if (is.null(existing)) {
        return(covariates)
      } else {
        return(unique(c(existing, covariates)))
      }
    },

    #' Classifies variables into mediators and outcomes based on LHS/RHS usage.
    classify_variables = function(lhs_vars, rhs_vars, regression_info) {
      unique_lhs <- unique(lhs_vars)
      unique_rhs <- unique(rhs_vars)
      mediators <- intersect(unique_lhs, unique_rhs)
      outcomes <- setdiff(unique_lhs, unique_rhs)

      output <- list(mediators = list(), outcomes = list())
      for (var in mediators) {
        output$mediators[[var]] <- regression_info[[var]]
      }
      for (var in outcomes) {
        output$outcomes[[var]] <- regression_info[[var]]
      }

      return(output)
    }
  )
)
