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

test_that("model_loader ignores comments and warns on unused lines", {
  vars <- mock_var_names_config$new()
  commented_txt <- "
    med =~ i1 + i2 + i3
    # This is a comment: out ~ i5
    out ~ treat + med + i4 # trailing comment"
  expect_no_warning(
    loader <- model_loader$new(commented_txt, vars, mock_cfa_model,
                               mock_reg_model)
  )
  expect_equal(loader$reg$txt, "out ~ treat + med + i4")

  intercept_txt <- "
    med =~ i1 + i2 + i3
    out ~ treat + med + i4
    i4 ~ 1"
  expect_warning(
    model_loader$new(intercept_txt, vars, mock_cfa_model, mock_reg_model),
    "not used by RAPSEM"
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

test_that("create_group_lavaan handles syntax correctly", {
  lavaan_txt <- paste(
    "f1 =~ i1 + i2",
    "med ~ treat",
    "out ~ treat + med + cov1",
    sep = "\n"
  )
  vars <- mock_var_names_config$new()
  model <- model_loader$new(lavaan_txt, vars)

  result_lines <- unlist(strsplit(model$group_lavaan, "\n"))
  expect_true(
    "f1 =~ i1 + i2" %in% result_lines
  )
  expect_true(
    "med ~ c(p_med_0, p_med_1)*med" %in% result_lines == FALSE
  )
  expect_true(
    "out ~ c(p_out_med_0, p_out_med_1)*med + c(p_out_cov1_0, p_out_cov1_1)*cov1"
    %in% result_lines
  )
})

test_that("residual covariances between indicators go to the CFA", {
  vars <- mock_var_names_config$new()
  txt <- "
    med =~ i1 + i2 + i3
    out ~ treat + med
    i1 ~~ i2"
  loader <- model_loader$new(txt, vars, mock_cfa_model, mock_reg_model)
  expect_true("i1 ~~ i2" %in% loader$cfa$txt)
  expect_false(any(grepl("~~", loader$reg$txt)))
})

test_that("latent covariances are not parsed as regressions", {
  vars <- mock_var_names_config$new()
  txt <- "
    med =~ i1 + i2 + i3
    out =~ i7 + i8 + i9
    out ~ treat + med
    med ~~ out"
  expect_warning(
    loader <- model_loader$new(txt, vars, mock_cfa_model, mock_reg_model),
    "not used by RAPSEM"
  )
  expect_equal(loader$reg$txt, "out ~ treat + med")
})

test_that("treatment by observed covariate interactions pass validation", {
  vars <- mock_var_names_config$new()
  txt <- "
    med =~ i1 + i2 + i3
    med ~ treat + cov2 + treat:cov2
    out ~ treat + med"
  expect_no_error(model_loader$new(txt, vars, mock_cfa_model, mock_reg_model))
})
