#' Empiricial Quantile Confidence Intervals
#'
#' @param x A `segmented_reg()` object.
#' @param reps Number of replicates.
#' @param conf.level The confidence level to use for the confidence interval.
#'   Must be strictly greater than 0 and less than 1.
#'   Defaults to 0.95 which corresponds to a 95 percent confidence interval.
#'
#' @return `x` but with empirical quantile confidence intervals instead of
#'   parametric confidence intervals (i.e., for APC and AAPC).
empirical_quantile_ci <- function(x, reps = 10000, conf.level = 0.95) {
  # 1. Estimate residuals from the original model.
  # 2. Compute the empirical CDF of residuals.
  # 3.	Generate uniform random values and obtain resampled residuals via inverse-CDF sampling.
  # 4.	Generate new outcome data by adding these resampled residuals to the model’s fitted values.
  # 5.	Refit the model to the resampled data and re-estimate AAPC.
  # 6.	Repeat many times and take quantiles of AAPC estimates for the CI.
  fit <- x$best_fit$model
  knots <- x$best_fit$metrics$knots
  n <- stats::nobs(fit)
  y_orig <- x$y
  x_orig <- x$x
  y_res <- stats::residuals(fit)
  y_est <- stats::fitted(fit)

  ecdf_res <- stats::ecdf(y_res)
  inv_ecdf_sample <- function(u, res) stats::quantile(res, probs = u, type = 1)

  result <- list()

  for (i in 1:reps) {
    u <- stats::runif(n)
    resampled_res <- inv_ecdf_sample(u, y_res)
    y_star <- exp(y_est + resampled_res)
    fit_star <- stats::lm(log(y_star) ~ lspline::lspline(x_orig, knots))
    coef_star <- stats::coef(fit_star)[-1]
    result_i <- 100 * (exp(coef_star) - 1)
    result[[i]] <- tibble::enframe(result_i)
  }

  emp_cis <-
    purrr::list_rbind(result) |>
    dplyr::summarise(
      nreps = n(),
      apc_ci_lwr = stats::quantile(value, probs = c((1 - conf.level) / 2)),
      apc_ci_upr = stats::quantile(value, probs = c(1 - (1 - conf.level) / 2)),
      .by = name
    )

  new_apc <-
    extract_best_apc(x) |>
    dplyr::select(-dplyr::contains('ci')) |>
    dplyr::bind_cols(emp_cis) |>
    dplyr::select(-name)

  x$best_fit$APC <- new_apc
  x$opts$conf.type <- "empirical quantile"
  x
}
