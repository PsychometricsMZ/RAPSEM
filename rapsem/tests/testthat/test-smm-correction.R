source(test_path("helpers/mock_helpers.R"))

test_that("create_correction_element delegates properly", {
  fs <- mock_factor_scores$new()
  fn <- fs$names
  sig <- fs$sigma
  preds <- cbind(mock_data, latent_vars)

  expect_equal(create_correction_element("cov1", preds, "med", preds,
                                         fn, sig),
               list(val = preds[, "cov1"] * preds[, "med"],
                    corr = - sig["cov1", "med"]))
  expect_equal(create_correction_element("cov1:cov2", preds, "med", preds,
                                         fn, sig),
               list(val = preds[, "cov1"] * preds[, "cov2"] * preds[, "med"],
                    corr = - preds[, "cov2"] * sig["cov1", "med"]))
  expect_equal(create_correction_element("med", preds, "cov1:cov2", preds,
                                         fn, sig),
               list(val = preds[, "cov1"] * preds[, "cov2"] * preds[, "med"],
                    corr = - preds[, "cov2"] * sig["cov1", "med"]))
  expect_equal(create_correction_element("cov1:cov2", preds, "med:cov2", preds,
                                         fn, sig),
               list(val = preds[, "cov1"] * preds[, "cov2"]
                    * preds[, "med"] * preds[, "cov2"],
                    corr = - preds[, "cov2"] * preds[, "cov2"]
                    * sig["cov1", "med"]))
})

test_that("handle_main_effect returns corrected element", {
  fs <- mock_factor_scores$new()
  fn <- fs$names
  sig <- fs$sigma
  preds <- cbind(mock_data, latent_vars)

  expect_equal(handle_main_effect("cov1", preds, "med", preds, fn, sig),
               list(val = preds[, "cov1"] * preds[, "med"],
                    corr =  - sig["cov1", "med"]))
  expect_equal(handle_main_effect("cov1", preds, "cov2", preds, fn, sig),
               list(val = preds[, "cov1"] * preds[, "cov2"],
                    corr = 0))
  expect_equal(handle_main_effect("cov2", preds, "treat", preds, fn, sig),
               list(val = preds[, "cov2"] * preds[, "treat"],
                    corr = 0))
})

test_that("extract_interaction_parts splits correctly", {
  fs <- mock_factor_scores$new()
  fn <- fs$names

  parts <- extract_interaction_parts("cov1:cov2", fn)
  expect_equal(parts$latent, "cov1")
  expect_equal(parts$observed, "cov2")
  expect_error(extract_interaction_parts("cov1:med", fn))
  parts <- extract_interaction_parts("cov2:treat", fn)
  expect_equal(parts$observed, c("cov2", "treat"))
  expect_error(extract_interaction_parts("cov1:cov2:med", fn))
})

test_that("handle_interaction returns corrected element", {
  fs <- mock_factor_scores$new()
  fn <- fs$names
  sig <- fs$sigma
  preds <- cbind(mock_data, latent_vars)

  expect_equal(handle_interaction("cov1:cov2", preds, "med", preds, fn, sig),
               list(val = preds[, "cov1"] * preds[, "cov2"] * preds[, "med"],
                    corr = - preds[, "cov2"] * sig["cov1", "med"]))
  expect_equal(handle_interaction("cov1:cov2:treat", preds, "med", preds,
                                  fn, sig),
               list(val = preds[, "cov1"] * preds[, "cov2"]
                    * preds[, "treat"] * preds[, "med"],
                    corr = - preds[, "cov2"] * preds[, "treat"]
                    * sig["cov1", "med"]))
  expect_equal(handle_interaction("cov1:cov2", preds, "cov2", preds, fn, sig),
               list(val = preds[, "cov1"] * preds[, "cov2"] * preds[, "cov2"],
                    corr = 0))
  expect_equal(handle_interaction("cov2:treat", preds, "cov2", preds, fn, sig),
               list(val = preds[, "cov2"] * preds[, "treat"] * preds[, "cov2"],
                    corr = 0))
})

test_that("handle_double_interaction returns corrected element", {
  fs <- mock_factor_scores$new()
  fn <- fs$names
  sig <- fs$sigma
  preds <- cbind(mock_data, latent_vars)

  expect_equal(handle_double_interaction("cov1:cov2", preds, "med:cov2", preds,
                                         fn, sig),
               list(val = preds[, "cov1"] * preds[, "cov2"]
                    * preds[, "med"] * preds[, "cov2"],
                    corr = - preds[, "cov2"] * preds[, "cov2"]
                    * sig["cov1", "med"]))
  expect_equal(handle_double_interaction("cov1:cov2", preds,
                                         "cov2:treat", preds,
                                         fn, sig),
               list(val = preds[, "cov1"] * preds[, "cov2"]
                    * preds[, "cov2"] * preds[, "treat"],
                    corr = 0))
  expect_equal(handle_double_interaction("treat:cov2", preds,
                                         "cov1:cov2", preds,
                                         fn, sig),
               list(val = preds[, "cov1"] * preds[, "cov2"]
                    * preds[, "cov2"] * preds[, "treat"],
                    corr = 0))
})
