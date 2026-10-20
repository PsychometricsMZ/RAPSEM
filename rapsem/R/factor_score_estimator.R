#' @title Factor Score Estimator Class
#' @description
#' Internal R6 class for computing factor scores from lavaan cfa output
#' according to the method of moments technique of Wall and Amemiya (2003).
#' @details
#' Takes a `cfa_model` object and a data sample as input, uses the lavaan
#' `cfa` function with ML estimation to fit the CFA model to the data
#' and computes the factor scores as well as the the error covariance matrix
#' based on extracted intercepts, loadings, and residual variances with
#' a method of  moments approach. McDonald's omega reliability coefficients
#' are also computed per latent factor.
#'
#' Missing indicator values: the CFA is fitted with full-information ML and
#' the scores are computed per missingness pattern from the observed
#' indicators. The Wall-Amemiya score equals the GLS score
#' (L' Psi^-1 L)^-1 L' Psi^-1 (z - mu) with error covariance
#' (L' Psi^-1 L)^-1; for a pattern with observed indicators O the same
#' formula is applied to L_O, Psi_OO and z_O, which gives a row-specific
#' error covariance (field `sigma_rows`).
#' @importFrom R6 R6Class
#' @importFrom lavaan cfa inspect coef
#' @keywords internal
factor_score_estimator <- R6Class("factor_score_estimator",
  public = list(
    #' @field names Character vector of latent factor names.
    names = NULL,
    #' @field vals Matrix of computed factor scores for each used observation.
    vals = NULL,
    #' @field sigma Covariance matrix of factor score estimation errors
    #'   (average over observations if they differ by missingness pattern).
    sigma = NULL,
    #' @field sigma_rows Array (observations x factors x factors) of
    #'   row-specific error covariances; NULL for complete indicator data.
    sigma_rows = NULL,
    #' @field rows Logical vector over the rows of the data sample: the rows
    #'   used for the scores and all subsequent regressions.
    rows = NULL,
    #' @field reliability Named numeric vector of McDonald's omega reliability
    #'   coefficients (one per latent factor)
    reliability = NULL,
    #' @field standardized Logical; TRUE if the scores and error covariances
    #'   are in the latent-standardized metric, FALSE if they are in the
    #'   metric of the anchor indicators (Wall & Amemiya).
    standardized = FALSE,

    #' @description
    #' Initializes a `factor_score_estimator` object and computes factor scores
    #' from a CFA model and data sample.
    #' @param cfa_model An instance of `cfa_model` defining the
    #'                  factor structure.
    #' @param data_sample A data frame containing observed indicator values.
    #' @param orthogonal Boolean which is True for orthogonal latent factors.
    #' @param keep_vars Observed (non-indicator) variables of the structural
    #'                  model; rows with missing values in them are not used.
    #' @param blockwise Boolean; if TRUE, the measurement model of each latent
    #'                  variable is fitted separately.
    #' @param standardize Boolean; if TRUE, the scores and error covariances
    #'                    are rescaled to the latent-standardized metric
    #'                    (latent variances 1). Default FALSE: metric of the
    #'                    anchor indicators, as in Wall & Amemiya.
    initialize = function(cfa_model, data_sample, orthogonal,
                          keep_vars = character(), blockwise = FALSE,
                          standardize = FALSE) {
      tryCatch({
        self$names <- cfa_model$fact_names
        self$standardized <- standardize
        results <- private$calculate_factor_scores(cfa_model, data_sample,
                                                   orthogonal, keep_vars,
                                                   blockwise, standardize)
        self$vals <- results$scores
        self$sigma <- results$sigma_ee
        self$sigma_rows <- results$sigma_rows
        self$rows <- results$rows
        self$reliability <- results$reliability
      }, error = function(e) {
        stop(paste("Error initializing factor score estimator: ",
                   conditionMessage(e)))
      })
    }
  ),

  private = list(
    #' Core logic for estimating factor scores from a CFA model and data sample.
    calculate_factor_scores = function(cfa_model, data, orthogonal,
                                       keep_vars = character(),
                                       blockwise = FALSE,
                                       standardize = FALSE) {
      free_indicators <- unlist(cfa_model$free_ind, use.names = FALSE)
      fixed_indicators <- unlist(cfa_model$fixed_ind, use.names = FALSE)
      all_indicators <- c(free_indicators, fixed_indicators)
      z_values <- data[, all_indicators]
      cfa_coef <- if (blockwise) {
        private$get_blockwise_cfa_coef(cfa_model, z_values)
      } else {
        private$get_cfa_coef(cfa_model, z_values, orthogonal)
      }
      beta0hat <- matrix(cfa_coef[paste0(free_indicators, "~1")],
                         nrow = length(free_indicators), ncol = 1)
      beta1hat <- private$build_loading_matrix(cfa_coef, cfa_model,
                                               free_indicators)
      psi_hat <- private$build_residual_covariance_matrix(cfa_coef,
                                                          all_indicators)
      num_facts <- length(cfa_model$fact_names)
      sigma_rows <- NULL
      if (!anyNA(z_values)) {
        b1_mod1 <- cbind(matrix(0, num_facts, length(free_indicators)),
                         diag(num_facts))
        b1_mod2 <- rbind(diag(length(free_indicators)), -1 * t(beta1hat))
        gamma <- (b1_mod1 %*% psi_hat %*% b1_mod2 %*%
                    solve(t(b1_mod2) %*% psi_hat %*% b1_mod2))
        sigma_ee <- (cbind(-gamma, diag(num_facts) + gamma %*% beta1hat)
                     %*% psi_hat %*% t(b1_mod1))
        num_obs <- nrow(data)
        fact_scores <- matrix(0, nrow = num_obs, ncol = num_facts)
        for (i in 1:num_obs) {
          centered <- t(as.matrix(z_values[i, ])) -
            rbind(beta0hat, matrix(0, num_facts, 1))
          fact_scores[i, ] <- (cbind(-gamma,
                                     diag(num_facts) + gamma %*% beta1hat)
                               %*% centered)
        }
      } else {
        patterned <- private$pattern_scores(
          as.matrix(z_values), rbind(beta1hat, diag(num_facts)),
          c(beta0hat, rep(0, num_facts)), psi_hat
        )
        fact_scores <- patterned$scores
        sigma_rows <- patterned$sigma_rows
      }
      rows <- stats::complete.cases(fact_scores)
      if (length(keep_vars) > 0) {
        rows <- rows & stats::complete.cases(data[, keep_vars, drop = FALSE])
      }
      if (!any(rows)) {
        stop("No observation has usable factor scores and complete data.")
      }
      fact_scores <- fact_scores[rows, , drop = FALSE]
      colnames(fact_scores) <- cfa_model$fact_names
      if (!is.null(sigma_rows)) {
        sigma_rows <- sigma_rows[rows, , , drop = FALSE]
        dimnames(sigma_rows) <- list(NULL, cfa_model$fact_names,
                                     cfa_model$fact_names)
        sigma_ee <- apply(sigma_rows, c(2, 3), mean)
      }
      colnames(sigma_ee) <- rownames(sigma_ee) <- cfa_model$fact_names
      if (standardize) {
        standardized <- private$standardize_latent(fact_scores, sigma_ee,
                                                   sigma_rows)
        fact_scores <- standardized$scores
        sigma_ee <- standardized$sigma_ee
        sigma_rows <- standardized$sigma_rows
      }
      reliability <- private$compute_reliability(cfa_coef, cfa_model,
                                                 beta1hat, free_indicators)
      list(
        scores      = fact_scores,
        sigma_ee    = sigma_ee,
        sigma_rows  = sigma_rows,
        rows        = rows,
        reliability = reliability
      )
    },

    #' Factor scores and error covariances per missingness pattern:
    #' GLS scores from the observed indicators of each row. Rows in which a
    #' factor has no observed indicator get NA scores. Residual variances
    #' fixed at 0 (see get_cfa_coef) are replaced by a tiny positive value,
    #' which gives the limit of the GLS score (that indicator determines the
    #' score of its factor) up to numerical precision.
    pattern_scores = function(z, loadings, intercepts, psi) {
      zero <- diag(psi) <= 0
      if (any(zero)) {
        diag(psi)[zero] <- 1e-10 * max(diag(psi))
      }
      num_facts <- ncol(loadings)
      scores <- matrix(NA_real_, nrow(z), num_facts)
      sigma_rows <- array(NA_real_, c(nrow(z), num_facts, num_facts))
      observed <- !is.na(z)
      pattern <- apply(observed, 1, function(o) paste(as.integer(o),
                                                      collapse = ""))
      for (p in unique(pattern)) {
        rows_p <- which(pattern == p)
        obs <- observed[rows_p[1], ]
        l_obs <- loadings[obs, , drop = FALSE]
        if (!any(obs) || qr(l_obs)$rank < num_facts) next
        psi_inv <- solve(psi[obs, obs, drop = FALSE])
        sigma_p <- solve(t(l_obs) %*% psi_inv %*% l_obs)
        weights <- sigma_p %*% t(l_obs) %*% psi_inv
        centered <- sweep(z[rows_p, obs, drop = FALSE], 2, intercepts[obs])
        scores[rows_p, ] <- centered %*% t(weights)
        for (r in rows_p) sigma_rows[r, , ] <- sigma_p
      }
      list(scores = scores, sigma_rows = sigma_rows)
    },

    #' Optional (standardize = TRUE): center the factor scores and scale them
    #' by the estimated latent standard deviation, sqrt(var(score) -
    #' sigma_ee), rescaling sigma_ee (and the row-specific sigma_rows)
    #' accordingly. Scores and error covariance stay in the same metric, and
    #' the corrected latent variances equal 1.
    standardize_latent = function(scores, sigma_ee, sigma_rows = NULL) {
      latent_var <- apply(scores, 2, stats::var) - diag(sigma_ee)
      if (any(!is.finite(latent_var) | latent_var <= 0)) {
        stop(paste0("Estimated latent variance is not positive for: ",
                    paste(colnames(scores)[latent_var <= 0], collapse = ", "),
                    ". The factor scores are too unreliable to be corrected."))
      }
      latent_sd <- sqrt(latent_var)
      centered <- sweep(scores, 2, colMeans(scores))
      scale_ee <- outer(latent_sd, latent_sd)
      if (!is.null(sigma_rows)) {
        for (i in seq_len(dim(sigma_rows)[1])) {
          sigma_rows[i, , ] <- sigma_rows[i, , ] / scale_ee
        }
      }
      list(
        scores     = sweep(centered, 2, latent_sd, "/"),
        sigma_ee   = sigma_ee / scale_ee,
        sigma_rows = sigma_rows
      )
    },

    #' Compute McDonald's omega per latent factor over all of its indicators
    #' (the fixed anchor indicator has loading 1):
    #' omega = (sum lambda)^2 psi / ((sum lambda)^2 psi + sum Theta), where
    #' psi is the latent variance and Theta the residual (co)variances of the
    #' factor's indicators.
    compute_reliability = function(cfa_coef, cfa_model, beta1hat,
                                   free_indicators) {
      reliability <- setNames(numeric(length(cfa_model$fact_names)),
                              cfa_model$fact_names)
      row_start <- 1
      for (j in seq_along(cfa_model$fact_names)) {
        fact        <- cfa_model$fact_names[j]
        indicators  <- cfa_model$free_ind[[fact]]
        num_ind     <- length(indicators)
        row_end     <- row_start + num_ind - 1

        loadings    <- c(1, beta1hat[row_start:row_end, j])
        all_ind     <- c(cfa_model$fixed_ind[[fact]], indicators)
        theta       <- private$build_residual_covariance_matrix(cfa_coef,
                                                                all_ind)
        psi         <- cfa_coef[paste0(fact, "~~", fact)]

        common      <- sum(loadings)^2 * psi
        reliability[fact] <- unname(common / (common + sum(theta)))

        row_start   <- row_end + 1
      }
      reliability
    },

    #' Fit cfa model, check that it converged and return coefficients.
    #' Negative variance estimates (Heywood cases) are handled following
    #' Kolenikov and Bollen (2012): a negative variance that is significantly
    #' below zero (Wald z < -1.96) indicates a misspecified model, and the
    #' solution is rejected (in the bootstrap, the replicate is skipped and
    #' counted). A non-significant negative residual variance is sampling
    #' error around a small positive value: it is fixed at its boundary 0
    #' and the model refitted (constrained ML); the indicator then determines
    #' the Wall-Amemiya score of its factor. The latent covariance matrix is
    #' not used by the scores, so it is not checked.
    get_cfa_coef = function(model, indicator_matrix, orthogonal) {
      syntax <- model$lavaan_syntax
      fit_cfa <- function(syntax, se = "none") {
        fit <- suppressWarnings(
          cfa(
            syntax,
            data = indicator_matrix,
            meanstructure = TRUE,
            estimator = "ML",
            missing = if (anyNA(indicator_matrix)) "ml" else "listwise",
            se = se,
            orthogonal = orthogonal,
            control = list(iter.max = 10000)
          )
        )
        if (!inspect(fit, "converged")) {
          stop("CFA model did not converge; cannot compute factor scores.")
        }
        fit
      }
      fit <- fit_cfa(syntax)
      indicators <- lavaan::lavNames(fit, "ov")
      # every refit fixes at least one more residual variance at 0
      for (iter in seq_len(length(indicators) + 1)) {
        pe <- lavaan::parameterEstimates(fit, se = FALSE, ci = FALSE)
        neg <- pe$op == "~~" & pe$lhs == pe$rhs & !is.na(pe$est) &
          pe$est < 0
        if (!any(neg)) break
        pe_se <- lavaan::parameterEstimates(fit_cfa(syntax, se = "standard"),
                                            ci = FALSE)
        z <- pe_se$est[neg] / pe_se$se[neg]
        clearly <- !is.finite(z) | z < stats::qnorm(0.025)
        if (any(clearly)) {
          stop("Improper CFA solution (Heywood case): variance of ",
               paste0(pe$lhs[neg][clearly], collapse = ", "),
               " significantly negative; cannot compute factor scores.")
        }
        fix <- pe$lhs[neg][pe$lhs[neg] %in% indicators]
        if (length(fix) == 0) break   # only latent variances, not used
        syntax <- private$fix_residual_variances(syntax, fix)
        fit <- fit_cfa(syntax)
      }
      theta <- lavaan::lavInspect(fit, "theta")
      if (min(eigen(theta, symmetric = TRUE, only.values = TRUE)$values) <
            -1e-8 * max(abs(diag(theta)))) {
        stop("Improper CFA solution: the residual covariance matrix is not ",
             "positive semi-definite; cannot compute factor scores.")
      }
      coef(fit)
    },

    #' Fix the residual variances of the given indicators at 0 in the lavaan
    #' syntax (replacing lines that free them explicitly).
    fix_residual_variances = function(syntax, indicators) {
      lines <- unlist(strsplit(syntax, "\n", fixed = TRUE))
      for (ind in indicators) {
        own <- grepl(paste0("^\\s*", ind, "\\s*~~\\s*", ind, "\\s*$"), lines)
        lines <- c(lines[!own], paste0(ind, " ~~ 0*", ind))
      }
      paste(lines, collapse = "\n")
    },

    #' Fit the measurement model of each latent variable separately and
    #' combine the coefficients. The Wall-Amemiya scores do not use the factor
    #' covariances, so they are not needed; misspecification in one block
    #' cannot affect the estimates of another.
    get_blockwise_cfa_coef = function(model, indicator_matrix) {
      lines <- trimws(unlist(strsplit(model$lavaan_syntax, "\n", fixed = TRUE)))
      lines <- lines[lines != ""]
      block_of <- function(line) {
        if (grepl("=~", line)) {
          return(trimws(sub("=~.*$", "", line)))
        }
        vars <- trimws(unlist(strsplit(gsub("~~", "+", line, fixed = TRUE),
                                       "+", fixed = TRUE)))
        owner <- vapply(vars, function(v) {
          hit <- Filter(function(f) v %in% c(model$fixed_ind[[f]],
                                             model$free_ind[[f]]),
                        model$fact_names)
          if (length(hit) == 1) hit else NA_character_
        }, character(1))
        if (anyNA(owner) || length(unique(owner)) > 1) {
          stop("Residual covariance '", line, "' links indicators of ",
               "different latent variables; use blockwise = FALSE.")
        }
        owner[[1]]
      }
      owners <- vapply(lines, block_of, character(1))
      coefs <- lapply(model$fact_names, function(fact) {
        inds <- c(model$fixed_ind[[fact]], model$free_ind[[fact]])
        if (length(inds) < 3) {
          stop("Latent variable '", fact, "' has fewer than 3 indicators and ",
               "cannot be identified on its own; use blockwise = FALSE.")
        }
        block <- list(lavaan_syntax = lines[owners == fact])
        private$get_cfa_coef(block, indicator_matrix[, inds, drop = FALSE],
                             orthogonal = FALSE)
      })
      unlist(coefs)
    },

    #' Extract factor loading matrix from lavaan fit
    build_loading_matrix = function(cfa_coef, cfa_model, free_indicators) {
      factor_loadings <- matrix(0, nrow = length(free_indicators),
                                ncol = length(cfa_model$fact_names))
      row_start <- 1
      for (j in seq_along(cfa_model$fact_names)) {
        fact <- cfa_model$fact_names[j]
        indicators <- cfa_model$free_ind[[fact]]
        loading_names <- paste0(fact, "=~", indicators)
        num_indicators <- length(indicators)
        row_end <- row_start + num_indicators - 1
        factor_loadings[row_start:row_end, j] <- cfa_coef[loading_names]
        row_start <- row_end + 1
      }
      factor_loadings
    },

    #' Extract residual covariance matrix from lavaan fit
    build_residual_covariance_matrix = function(cfa_coef, indicators) {
      cov_labels <- outer(indicators, indicators,
        FUN = Vectorize(function(i, j) {
          if (i <= j) paste0(i, "~~", j) else paste0(j, "~~", i)
        })
      )
      cov_values <- matrix(0, nrow = length(indicators),
                           ncol = length(indicators))
      for (i in seq_along(indicators)) {
        for (j in seq_along(indicators)) {
          label <- cov_labels[i, j]
          if (label %in% names(cfa_coef)) {
            cov_values[i, j] <- cfa_coef[label]
          }
        }
      }
      cov_values <- (cov_values + t(cov_values)) / 2
      colnames(cov_values) <- rownames(cov_values) <- indicators
      cov_values
    }
  )
)