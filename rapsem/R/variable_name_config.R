#' @title Variable Name Configuration Class
#' @description
#' Internal R6 class for storing variable name configuration.
#' Holds treatment names and all variable names for internal reference.
#' @importFrom R6 R6Class
#' @keywords internal
var_name_config <- R6Class("var_name_config",
  public = list(
    #' @field intercept String for intercept name; defaults to "intercept".
    intercept = NULL,
    #' @field treatment Character vector of treatment variable names.
    treatment = NULL,
    #' @field all Character vector of all variable names in the input data.
    all = NULL,

    #' @description
    #' Initialize the var_name_config object.
    #' @param treatment_names Character vector of treatment variable names.
    #' @param all_var_names Character vector of all variable names.
    #' @param intercept_name String for intercept name; defaults to "intercept".
    initialize = function(treatment_names, all_var_names,
                          intercept_name = "intercept") {
      self$intercept <- intercept_name
      private$validate_treatment_names(treatment_names, all_var_names)
      self$treatment <- treatment_names
      self$all <- all_var_names
    }
  ),

  private = list(
    #' Validates the treatment names against the input data column names.
    #' Ensures that at least one treatment is provided and that no
    #' treatment is the same as the intercept name.
    validate_treatment_names = function(treatment_names, all_var_names) {
      if (length(treatment_names) == 0) {
        stop("At least one treatment variable must be specified.")
      }
      for (treat_name in treatment_names) {
        if (!(treat_name %in% all_var_names)) {
          stop(paste("treatment", treat_name, "must be in column names."))
        }
      }
      if (self$intercept %in% treatment_names) {
        stop(paste(self$intercept,
                   "cannot be used as a treatment variable name."))
      }
    }
  )
)
