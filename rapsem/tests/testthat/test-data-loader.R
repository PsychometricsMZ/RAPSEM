test_that("data_loader initializes correctly with valid data.frame input", {
  test_data <- data.frame(treat = c(0, 1, 0, 1, 1), cov = seq(10, 50, by = 10))
  treatment_names <- "treat"

  loader <- data_loader$new(test_data, treatment_names)

  expect_s3_class(loader, "data_loader")
  expect_equal(loader$original_data[, c("treat", "cov")], test_data)
  expect_s3_class(loader$name_config, "var_name_config")
  expect_equal(loader$name_config$treatment, treatment_names)
  expect_equal(loader$name_config$all, colnames(test_data))
  intercept_name <- loader$name_config$intercept
  expect_true(intercept_name %in% colnames(loader$original_data))
  expect_equal(unique(loader$original_data[[intercept_name]]), 1)
})

test_that("data_loader handles multiple treatment names correctly", {
  df <- data.frame(treat1 = c(0, 1, 1, 0, 1),
                   treat2 = c(1, 0, 0, 1, 1),
                   y = 11:15)
  treatment_names <- c("treat1", "treat2")
  loader <- data_loader$new(df, treatment_names)

  expect_equal(loader$name_config$treatment, treatment_names)
  expect_true(all(treatment_names %in% colnames(loader$original_data)))
})

test_that("data_loader accepts numeric matrix and converts to data.frame", {
  m <- matrix(c(0, 1, 1, 2, 3, 4, 5, 6, 7), nrow = 3)
  colnames(m) <- c("treat", "cov", "missing")
  m[, "treat"] <- c(0, 1, 0)
  loader <- data_loader$new(m, "treat")
  expect_s3_class(loader$original_data, "data.frame")
  expect_equal(ncol(loader$original_data), 4)  # includes intercept
  expect_equal(colnames(loader$original_data), c("treat", "cov",
                                                 "missing", "intercept"))
  expect_equal(unique(loader$original_data$intercept), 1)
})

test_that("data_loader handles empty dataset (zero rows)", {
  empty_df <- data.frame(x = numeric(0), y = numeric(0))
  expect_error(data_loader$new(empty_df, "x"),
               "input_data must have at least one row.")
})

test_that("data_loader throws error if data has non-numeric columns", {
  df <- data.frame(x = 1:5, y = letters[1:5])
  expect_error(data_loader$new(df, "x"),
               "All columns in input_data must be numeric")
})

test_that("data_loader throws error if matrix is not numeric", {
  m <- matrix(as.character(1:9), nrow = 3)
  colnames(m) <- c("treat", "cov", "missing")
  expect_error(data_loader$new(m, "treat"),
               "All columns in input_data must be numeric")
})

test_that("data_loader throws error for invalid input type", {
  expect_error(data_loader$new(1:5, "treat"),
               "input_data must be a matrix or data.frame.")
})

test_that("data_loader throws error if input has no column names", {
  no_colnames_matrix <- matrix(1:9, nrow = 3)
  colnames(no_colnames_matrix) <- NULL
  expect_error(data_loader$new(no_colnames_matrix, "treat"),
               "input_data must have column names.")
})

test_that("data_loader throws error for duplicated column names", {
  dup_data <- matrix(1:9, nrow = 3)
  colnames(dup_data) <- c("treat", "cov", "treat")
  expect_error(data_loader$new(dup_data, "treat"),
               "input_data must not have duplicate column names.")
})

test_that("data_loader overwrites intercept column with warning", {
  df <- data.frame(treat = rep(c(0, 1), 5), intercept = rnorm(5))
  expect_warning(loader <- data_loader$new(df, "treat"),
                 "already exists and will be overwritten")
  expect_equal(unique(loader$original_data$intercept), 1)
  expect_true(all(loader$original_data$intercept == 1))
})

test_that("data_loader throws error if treatment not in columns", {
  df <- data.frame(x = 1:5, y = 1:5)
  expect_error(data_loader$new(df, "z"),
               "treatment z must be in column names")
})

test_that("data_loader throws error if intercept is used as treatment", {
  df <- data.frame(intercept = rep(c(0, 1), 5), x = rep(c(0, 1), 5))
  expect_error(data_loader$new(df, "intercept"),
               "intercept cannot be used as a treatment")
})

test_that("data_loader throws error if any treatment name is missing", {
  df <- data.frame(treat1 = rep(c(0, 1), 5), y = 6:10)
  expect_error(
    data_loader$new(df, c("treat1", "treat2")),
    "treatment treat2 must be in column names"
  )
})

test_that("data_loader throws error if intercept is among treatment names", {
  df <- data.frame(x = 1:5, y = 6:10, intercept = 11:15)
  expect_error(
    data_loader$new(df, c("x", "intercept")),
    "intercept cannot be used as a treatment"
  )
})

test_that("data_loader throws error if no treatment names provided", {
  df <- data.frame(x = 1:5, y = 6:10)
  expect_error(
    data_loader$new(df, character(0)),
    "At least one treatment variable must be specified."
  )
})

test_that("data_loader throws error when treatment columns are not binary", {
  df <- data.frame(
    treat1 = c(0, 1, 0, 2, 1),
    treat2 = c(0, 1, 0, 1, 3),
    cov = rnorm(5)
  )
  expect_error(
    data_loader$new(df, c("treat1", "treat2")),
    paste("Treatment columns must contain only 0 or 1.",
          "Invalid columns: treat1, treat2")
  )
})

test_that("bootstrap generates correct number of bootstrap samples", {
  set.seed(123)
  df <- data.frame(x = rep(0:1, 5), y = rnorm(10))
  loader <- data_loader$new(df, "x")
  loader$bootstrap(5)
  expect_length(loader$boot_samples, 5)
  expect_true(all(sapply(loader$boot_samples, nrow) == nrow(df)))
  expect_true(all(sapply(loader$boot_samples,
                         function(df_) {
                           all(colnames(df_) == colnames(loader$original_data))
                         })))
})

test_that("bootstrap works correctly with n_bootstrap = 1", {
  set.seed(123)
  df <- data.frame(x = rep(0:1, 5), y = rnorm(10))
  loader <- data_loader$new(df, "x")
  loader$bootstrap(1)
  expect_length(loader$boot_samples, 1)
  expect_equal(nrow(loader$boot_samples[[1]]), nrow(df))
})

test_that("bootstrap works correctly when n_bootstrap > n observations", {
  set.seed(123)
  df <- data.frame(x = rep(0:1, 5), y = rnorm(5))
  loader <- data_loader$new(df, "x")
  loader$bootstrap(10)
  expect_length(loader$boot_samples, 10)
  expect_true(all(sapply(loader$boot_samples, function(sample) {
    nrow(sample) == nrow(df)
  })))
})


test_that("bootstrap with n = 0 returns empty list", {
  df <- data.frame(x = rep(c(0, 1), 5), y = rnorm(5))
  loader <- data_loader$new(df, "x")
  loader$bootstrap(0)
  expect_equal(loader$boot_samples, list())
})
