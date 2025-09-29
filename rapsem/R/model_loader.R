#' @title Model Loader Class
#' @description
#' Internal R6 class to parse and validate lavaan model text,
#' and organize it into CFA and regression components.
#' @details
#' The model loader takes lavaan model syntax, validates it
#' against provided variable names, and separates it into
#' confirmatory factor analysis (CFA) and regression components,
#' each of which are stored in separate R6 classes
#' (`cfa_model`, `reg_model`).
#' @importFrom R6 R6Class
#' @importFrom lavaan lavaanify lavNames
#' @keywords internal
model_loader <- R6Class("model_loader",
  public = list(
    #' @field cfa An instance of the `cfa_model` class
    #' representing CFA components.
    cfa = NULL,
    #' @field reg An instance of the `reg_model` class
    #' representing regression components.
    reg = NULL,

    #' @description
    #' Initialize the model_loader object.
    #' @param lavaan_text A string containing lavaan model syntax.
    #' @param var_names An instance of `var_names_config`,
    #'                  containing the variable names.
    #' @param cfa_class A class generator for CFA models,
    #'                  allowing dependency injection.
    #' @param reg_class A class generator for regression models,
    #'                  allowing dependency injection.
    initialize = function(lavaan_text, var_names,
                          cfa_class = cfa_model,
                          reg_class = reg_model) {
      tryCatch({
        lavaan_text <- private$validate_lavaan_txt(lavaan_text, var_names$all)
        split_lavaan_txt <- private$split_lavaan_txt(lavaan_text)
        cfa_txt <- private$validate_lavaan_component(split_lavaan_txt$cfa,
                                                     "cfa")
        reg_txt <- private$validate_lavaan_component(split_lavaan_txt$reg,
                                                     "reg")
        self$cfa <- cfa_class$new(cfa_txt)
        self$reg <- reg_class$new(reg_txt, var_names$intercept)
      }, error = function(e) {
        stop("Error initializing model_loader: ", e$message)
      })
    }
  ),

  private = list(
    #' Validates the lavaan model text syntax and checks for variables
    #' not present in the data.
    validate_lavaan_txt = function(model_text, all_var_names) {
      if (is.character(model_text) && length(model_text) > 1) {
        model_text <- paste(model_text, collapse = "\n")
      }
      if (!requireNamespace("lavaan", quietly = TRUE)) {
        stop("The 'lavaan' package is required but not installed.")
      }
      if (!is.character(model_text) || length(model_text) != 1) {
        stop("lavaan_text must be a single character string.")
      }
      model_parsed <- tryCatch(
        lavaanify(model_text),
        error = function(e) stop("Invalid lavaan model syntax: ", e$message)
      )
      model_vars <- lavNames(model_parsed, type = "ov")
      missing_vars <- setdiff(model_vars, all_var_names)
      if (length(missing_vars) > 0) {
        stop("Variables in model not found in data: ",
             paste(missing_vars, collapse = ", "))
      }
      return(model_text)
    },

    #' Separates lavaan syntax lines into CFA and regression components.
    split_lavaan_txt = function(model_text) {
      lines <- private$split_txt_into_lines(model_text)
      cfa_lines <- c()
      regression_lines <- c()
      for (line in lines) {
        if (grepl("=~", line)) {
          cfa_lines <- c(cfa_lines, line)
        } else if (grepl("~", line)) {
          regression_lines <- c(regression_lines, line)
        } else {
          warning(paste("Line", line, "not considered."))
        }
      }
      return(list(cfa = cfa_lines, reg = regression_lines))
    },

    #' Trims and splits the model text into individual non-empty lines.
    split_txt_into_lines = function(text) {
      if (!is.character(text) || length(text) != 1) {
        stop("Input must be a single character string.")
      }
      lines <- unlist(strsplit(text, "\n"))
      lines <- trimws(lines)
      return(lines[lines != ""])
    },

    #' Checks that a lavaan submodel component is present.
    validate_lavaan_component = function(component_txt, type) {
      if (length(component_txt) == 0) {
        stop(paste(type, "component missing in lavaan model."))
      }
      return(component_txt)
    }
  )
)
