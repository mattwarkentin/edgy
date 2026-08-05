#' Segmented Regression using Linear Splines
#'
#' `segmented_reg()` fits a linear spline model for the outcome (e.g., rates)
#'   against a time variables (e.g., years). Grid search is performed to find
#'   the best model (i.e., optimal knot locations) according to the `metric`
#'   criterion. Information about the best model can be extracted from the
#'   returned object using the `extract_best_*()` set of functions.
#'
#' @param formula Model formula (e.g., `year ~ rate`).
#' @param data Data frame or [`tibble::tibble`] with variables in `formula`.
#' @param events Vector with the number of events (e.g., deaths) that correspond
#'   to the rate data in the `formula.` If provided, weighted least squares is
#'   used instead of ordinary least squares for fitting regression models.
#' @param opts A `knots_opt()` object providing knot location options.
#' @param metric Metric to use for model selection. One of `"bic"`, `"bic3"`,
#'   `"aic"`, or `"aicc".` Default is `"bic3"`.
#' @param conf_level The confidence level to use for the confidence interval.
#'   Must be strictly greater than 0 and less than 1.
#'   Defaults to 0.95 which corresponds to a 95 percent confidence interval.
#' @param progress Whether to show progress bars. By default, progress bars are
#'   enabled in interactive sessions (i.e., if `rlang::is_interactive()`
#'   returns `TRUE`).
#' @param knots A numeric vector of knot locations or a `list` of numeric
#'   vectors specifying sets of knot locations. If provided, this overrides the
#'   automatic knot selection performed using `knot_opts()` and `data`.
#' @param x A `segmented_reg()` object.
#' @param ... Not currently used.
#'
#' @details
#' [lspline::lspline()] is used to compute the basis for a piecewise linear
#'   spline to estimate coefficients in the segmented regression model.
#'
#'   If [mirai::daemons()] has been used to set persistent background processes,
#'   this function will fit segmented regression models in parallel using all
#'   available processes.
#'
#' @return A named-list of class `"edgy_segmented_reg"`.
#'
#' @md
#'
#' @importFrom rlang :=
#' @import carrier
#' @import mirai
#'
#' @examples
#' df <- read.delim(system.file("example.txt", package = "edgy"), header = FALSE)
#' res <- segmented_reg(V3 ~ V2, data = df)
#'
#' @md
#' @export
segmented_reg <- function(
  formula,
  data,
  events = NULL,
  opts = knot_opts(),
  metric = 'bic3',
  conf_level = 0.95,
  progress = rlang::is_interactive(),
  knots,
  ...
) {
  x_var <- rlang::f_rhs(formula)
  x_vals <- rlang::inject(`$`(data, !!x_var))
  y_var <- rlang::f_lhs(formula)
  y_vals <- rlang::inject(`$`(data, !!y_var))

  check_positive_rates(y_vals)

  if (rlang::is_null(events)) {
    wts <- rep_len(1L, length(y_vals))
  } else {
    check_nonnegative_int(events)
    wts <- events
  }

  metric_choices <- c('bic', 'bic3', 'aic', 'aicc')
  metric <- rlang::arg_match(metric, metric_choices)

  opts$fun <- log
  opts$inv_fun <- exp
  opts$conf.type <- "parametric"
  opts$metrics <- metric_choices

  no_knot_form <- rlang::inject(opts$fun(!!y_var) ~ !!x_var)

  no_knot_model <- stats::lm(formula = no_knot_form, data = data, weights = wts)

  no_knot_model <- structure(
    no_knot_model,
    class = c("edgy_spline_fit", class(no_knot_model))
  )

  no_knot_data <-
    tibble::enframe(
      list(no_knot_model),
      name = NULL,
      value = 'model'
    ) |>
    dplyr::mutate(nknots = 0)

  if (rlang::is_missing(knots)) {
    knot_outputs <- make_knot_sets(x_vals, opts)
    opts <- knot_outputs$opts
    knot_sets <- knot_outputs$knot_sets
  } else {
    if (rlang::is_bare_numeric(knots)) {
      knots <- list(knots)
    }
    check_list_of_numeric_vectors(knots)
    knot_sets <- knots
    no_knot_data <- NULL
  }

  res <-
    knot_sets |>
    tibble::enframe(name = NULL, value = "knots") |>
    dplyr::mutate(
      model = purrr::map(
        .x = knots,
        .f = purrr::in_parallel(
          \(k) {
            segmented_reg_fit(y_vals, x_vals, wts, k, opts)
          },
          segmented_reg_fit = segmented_reg_fit,
          y_vals = y_vals,
          x_vals = x_vals,
          wts = wts,
          opts = opts
        ),
        .progress = ifelse(progress, "Fitting models", FALSE)
      ),
      nknots = purrr::map_int(knots, length)
    ) |>
    dplyr::bind_rows(no_knot_data) |>
    dplyr::arrange(nknots) |>
    dplyr::mutate(
      id = 1:dplyr::n(),
      preds = purrr::map(model, \(m) m$fitted.values),
      k = purrr::map_int(model, \(m) base::length(m$coefficients) + 1),
      L = purrr::map_dbl(model, \(m) as.numeric(stats::logLik(m))),
      N = purrr::map_int(model, stats::nobs),
      df = get_deg_free(N, nknots),
      sse = purrr::map_dbl(model, \(m) sum(m$residuals^2)),
      aic = -2 * L + 2 * k,
      aicc = aic + 2 * (k * (k + 1)) / (N - k - 1),
      bic = log(sse / N) + ((2 * k + 2) / N) * log(N),
      bic3 = log(sse / N) + ((3 * k + 2) / N) * log(N)
    ) |>
    dplyr::select(id, dplyr::everything())

  best_fit <- dplyr::slice_min(res, !!rlang::sym(metric), with_ties = FALSE)

  data_with_preds <- dplyr::mutate(
    .data = data,
    est = opts$inv_fun(best_fit$preds[[1]])
  )

  apc_data <- estimate_apc(
    x = x_vals,
    knots = best_fit$knots[[1]],
    model = best_fit$model[[1]],
    deg_free = best_fit$df,
    conf_level = conf_level,
    opts = opts
  )

  apc <-
    apc_data |>
    dplyr::select(
      segment,
      period_start,
      period_end,
      apc,
      apc_ci_lwr,
      apc_ci_upr,
      apc_pval
    )

  aapc <- estimate_aapc(
    apc_data,
    deg_free = best_fit$df,
    conf_level = conf_level,
    opts = opts
  )

  ret <- list(
    formula = formula,
    y_var = y_var,
    x_var = x_var,
    y_vals = y_vals,
    x_vals = x_vals,
    data = data_with_preds,
    opts = opts,
    knot_sets = knot_sets,
    fits = res,
    best_fit = list(
      fit = best_fit,
      model = best_fit$model[[1]],
      criterion = metric,
      conf_level = conf_level,
      deg_free = best_fit$df,
      knots = best_fit$knots[[1]],
      metrics = list(
        nknots = best_fit$nknots,
        nobs = best_fit$N,
        logLik = best_fit$L,
        AIC = best_fit$aic,
        AICc = best_fit$aicc,
        BIC = best_fit$bic,
        BIC3 = best_fit$bic3
      ),
      APC = apc,
      AAPC = aapc
    )
  )
  structure(
    ret,
    class = c("edgy_segmented_reg", class(ret))
  )
}

segmented_reg_fit <- function(y, x, w, k, opts) {
  fit <- stats::lm(
    formula = rlang::inject(opts$fun(!!y) ~ lspline::lspline(!!x, k)),
    weights = w
  )
  structure(fit, class = c("edgy_spline_fit", class(fit)))
}

# https://surveillance.cancer.gov/help/joinpoint/statistical-notes/statistics-related-to-the-k-joinpoint-model/degrees-of-freedom
get_deg_free <- function(n, k) {
  (n - k) - (2 * (k + 1))
}

#' @rdname segmented_reg
#' @export
print.edgy_segmented_reg <- function(x, ...) {
  nmodels <- format(nrow(x$fits), big.mark = ',')
  mets <- x$best_fit$metrics
  locs <- glue::glue_collapse(x$best_fit$knots, ", ")
  cli::cat_rule(glue::glue('Segmented Regression ({nmodels} models evaluated)'))
  cli::cli_alert(glue::glue('Number of obs: {mets$nobs}'))
  cli::cli_alert(glue::glue('Number of knots: {mets$nknots}'))
  cli::cli_alert(glue::glue('Knot locations: {locs}'))
  invisible(x)
}
