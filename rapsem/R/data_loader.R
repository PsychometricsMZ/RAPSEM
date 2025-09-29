#' @title Data Loader Class
#' @description
#' Internal R6 class for loading and bootstrapping input data.
#' @details
#' Validates and stores input data, manages variable names,
#' and generates bootstrap samples.
#' @importFrom R6 R6Class
#' @keywords internal
data_loader <- R6Class("data_loader",
  public = list(
    #' @field original_data A matrix or data.frame of the input data.
    original_data = NULL,
    #' @field name_config An instance of the `var_name_config` class.
    name_config = NULL,
    #' @field boot_samples A list of bootstrap samples
    #'                     generated from `original_data`.
    boot_samples = NULL,

    #' @description
    #' Initialize the data_loader object.
    #' @param input_data A matrix or data.frame with named columns
    #'                   containing the input data.
    #' @param treatment_names A character vector of treatment variable names.
    initialize = function(input_data, treatment_names) {
      input_data <- private$validate_input_data(input_data)
      self$name_config <- var_name_config$new(treatment_names,
                                              colnames(input_data))
      private$validate_treatment_vars(input_data, treatment_names)
      intercept_col <- self$name_config$intercept
      if (intercept_col %in% colnames(input_data)) {
        warning(paste("Column", intercept_col, "already exists and will be",
                      "overwritten."))
      }
      input_data[, intercept_col] <- 1
      self$original_data <- input_data
    },

    #' @description
    #' Generate bootstrap samples from the original data.
    #' @param n_bootstrap Integer, number of bootstrap samples to generate.
    bootstrap = function(n_bootstrap) {
      if (!is.numeric(n_bootstrap) || length(n_bootstrap) != 1
          || n_bootstrap < 0 || !is.finite(n_bootstrap)) {
        stop("n_bootstrap must be a non-negative integer.")
      }
      self$boot_samples <- private$generate_bootstrap_samples(n_bootstrap)
    }
  ),

  private = list(
    #' Validates the format and column names of the input data.
    validate_input_data = function(input_data) {
      if (!inherits(input_data, c("matrix", "data.frame"))) {
        stop("input_data must be a matrix or data.frame.")
      }
      if (nrow(input_data) == 0) {
        stop("input_data must have at least one row.")
      }
      if (is.null(colnames(input_data))) {
        stop("input_data must have column names.")
      }
      if (anyDuplicated(colnames(input_data))) {
        stop("input_data must not have duplicate column names.")
      }
      if (inherits(input_data, "matrix")) {
        input_data <- as.data.frame(input_data)
      }
      if (!all(sapply(input_data, is.numeric))) {
        stop("All columns in input_data must be numeric.")
      }
      return(input_data)
    },

    #' Validates that treatment variables have binary format.
    validate_treatment_vars = function(input_data, treatment_names) {
      non_binary <- sapply(input_data[treatment_names], function(col) {
        !all(col %in% c(0, 1))
      })
      if (any(non_binary)) {
        stop(paste("Treatment columns must contain only 0 or 1.",
                   "Invalid columns:",
                   paste(names(non_binary)[non_binary], collapse = ", ")))
      }
    },

    # Generates bootstrap samples
    generate_bootstrap_samples = function(n_bootstrap) {
      if (n_bootstrap == 0) {
        return(list())
      }
      num_obs <- nrow(self$original_data)
      boot_samples <- vector("list", n_bootstrap)
      for (i in seq_len(n_bootstrap)) {
        sample_indices <- sample(seq_len(num_obs), size = num_obs,
                                 replace = TRUE)
        sample <- self$original_data[sample_indices, , drop = FALSE]
        colnames(sample) <- colnames(self$original_data)
        boot_samples[[i]] <- sample
      }
      return(boot_samples)
    }
  )
)
