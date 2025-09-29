source(test_path("helpers/mock_helpers.R"))

test_that("model_loader initializes correctly with valid input", {
  lavaan_txt <- "med =~ i1 + i2 + i3
                 out ~ treat + med + i4"
  vars <- mock_var_names_config$new()
  mdl <- model_loader$new(lavaan_txt, vars, mock_cfa_model, mock_reg_model)

  expect_true(inherits(mdl, "model_loader"))
  expect_true(inherits(mdl$cfa, "cfa_model"))
  expect_equal(trimws(mdl$cfa$txt), "med =~ i1 + i2 + i3")
  expect_true(inherits(mdl$reg, "reg_model"))
  expect_equal(trimws(mdl$reg$txt), "out ~ treat + med + i4")
})

test_that("model_loader handles multiple CFA and regression lines", {
  lavaan_txt <- "cov1 =~ i1 + i2 + i3
                 med =~ i4 + i5 + i6
                 out =~ i7 + i8 + i9
                 med ~ treat + cov2
                 out ~ treat + med + cov1"
  vars <- mock_var_names_config$new()
  mdl <- model_loader$new(lavaan_txt, vars, mock_cfa_model, mock_reg_model)

  expect_length(mdl$cfa$txt, 3)
  expect_length(mdl$reg$txt, 2)
})

test_that("model_loader warns on unrecognized syntax line", {
  lavaan_txt <- "
    med =~ i1 + i2 + i3
    # This is a comment
    out ~ treat + med + i4"
  vars <- mock_var_names_config$new()

  expect_warning(
    model_loader$new(lavaan_txt, vars, mock_cfa_model, mock_reg_model),
    "not considered"
  )
})

test_that("model_loader throws error with invalid lavaan syntax", {
  bad_txt <- "med =~ i1 + i2 + i3
              out ~ treat + med + i4
              out !~ i5"
  vars <- mock_var_names_config$new()

  expect_error(model_loader$new(bad_txt, vars),
               "Invalid lavaan model syntax")
})

test_that("model_loader throws error if variables are not in var_names", {
  txt <- "med =~ i1 + i2 + i3
          out ~ treat + med + x4"
  vars <- mock_var_names_config$new()

  expect_error(model_loader$new(txt, vars),
               "Variables in model not found in data")
})

test_that("model_loader throws error if CFA component is missing", {
  txt <- "out ~ treat + i4"
  vars <- mock_var_names_config$new()

  expect_error(model_loader$new(txt, vars),
               "cfa component missing in lavaan model")
})

test_that("model_loader throws error if regression component is missing", {
  txt <- "med =~ i1 + i2 + i3"
  vars <- mock_var_names_config$new()

  expect_error(model_loader$new(txt, vars),
               "reg component missing in lavaan model")
})
