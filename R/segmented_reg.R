#' Segmented Regression using Linear Splines
#'
#' `segmented_reg()` fits a linear spline model for the outcome (e.g., rates)
#'   against a time variables (e.g., years). Grid search is performed to find
#'   the best model (i.e., optimal knot locations) according to the `metric`
#'   criterion. Information about the best model can be extracted from the
#'   returned object using the `extract_best_*()` set of functions.
#'
#' @param formula Model formula.
#' @param data Data frame or `tibble`.
#' @param opts A `knots_opt()` object providing knot location options.
#' @param metric Metric to use for model selection. One of `"bic"`, `"bic3"`,
#'   `"aic"`, or `"aicc".` Default is `"bic3"`.
#' @param conf_level The confidence level to use for the confidence interval.
#'   Must be strictly greater than 0 and less than 1.
#'   Defaults to 0.95 which corresponds to a 95 percent confidence interval.
#' @param x A `segmented_reg()` object.
#' @param ... Not currently used.
#'
#' @details
#' `lspline::lspline()` is used to compute the basis for a piecewise linear
#'   spline to estimate coefficients in the segmented regression model.
#'
#' @return A named-list of class `"edgy_segmented_reg"`.
#'
#' @md
#'
#' @importFrom rlang :=
#'
#' @examples
#' \dontrun{
#' fit <- segmented_reg(rate ~ year, data = df)
#'
#' extract_best_model(fit)
#' extract_best_metrics(fit)
#' extract_best_predictions(fit)
#' extract_best_apc(fit)
#' extract_best_aapc(fit)
#' }
#'
#' @md
#' @export
segmented_reg <- function(
  formula,
  data,
  opts = knot_opts(),
  metric = 'bic3',
  conf_level = 0.95,
  ...
) {
  x_var <- rlang::f_rhs(formula)
  x_vals <- rlang::inject(`$`(data, !!x_var))
  y_var <- rlang::f_lhs(formula)
  y_vals <- rlang::inject(`$`(data, !!y_var))

  metric <- rlang::arg_match(metric, c('bic', 'bic3', 'aic', 'aicc'))

  opts$fun <- log
  opts$inv_fun <- exp
  opts$conf.type <- "parametric"

  knot_outputs <- make_knot_sets(x_vals, opts)
  opts <- knot_outputs$opts
  knot_sets <- knot_outputs$knot_sets

  no_knot_form <- rlang::inject(opts$fun(!!y_var) ~ !!x_var)

  no_knot_model <- stats::lm(formula = no_knot_form, data = data)

  no_knot_data <-
    tibble::enframe(
      list(no_knot_model),
      name = NULL,
      value = 'model'
    ) |>
    dplyr::mutate(nknots = 0)

  res <-
    knot_sets |>
    tibble::enframe(name = NULL, value = 'knots') |>
    dplyr::mutate(
      model = purrr::map(knots, \(k) {
        segmented_reg_fit(y_var, x_var, k, data, opts)
      }),
      nknots = purrr::map_int(knots, length)
    ) |>
    dplyr::bind_rows(no_knot_data) |>
    dplyr::mutate(
      preds = purrr::map(model, \(m) {
        suppressWarnings(stats::predict(m, interval = 'prediction'))
      }),
      preds = purrr::map(preds, tibble::as_tibble),
      k = purrr::map_int(model, \(m) base::length(m$coefficients) + 1),
      L = purrr::map_dbl(model, \(m) as.numeric(stats::logLik(m))),
      N = purrr::map_int(model, stats::nobs),
      sse = purrr::map_dbl(model, \(m) sum(m$residuals^2)),
      aic = -2 * L + 2 * k,
      aicc = aic + 2 * (k * (k + 1)) / (N - k - 1),
      bic = log(sse / N) + ((2 * k + 2) / N) * log(N),
      bic3 = log(sse / N) + ((3 * k + 2) / N) * log(N)
    ) |>
    dplyr::arrange(nknots)

  best_fit <- dplyr::slice_min(res, !!rlang::sym(metric), with_ties = FALSE)

  data_with_preds <-
    dplyr::bind_cols(
      data,
      best_fit |>
        dplyr::select(preds) |>
        tidyr::unnest(preds)
    ) |>
    dplyr::mutate(dplyr::across(c(fit, lwr, upr), opts$inv_fun)) |>
    dplyr::rename(
      est = fit
    )

  # https://surveillance.cancer.gov/help/joinpoint/statistical-notes/statistics-related-to-the-k-joinpoint-model/degrees-of-freedom
  deg_free <- (best_fit$N - best_fit$nknots) - (2 * (best_fit$nknots + 1))

  apc_data <- estimate_apc(
    x = x_vals,
    knots = best_fit$knots[[1]],
    model = best_fit$model[[1]],
    deg_free = deg_free,
    conf_level = conf_level,
    opts = opts
  )

  apc <-
    apc_data |>
    dplyr::select(
      period_start,
      period_end,
      apc,
      apc_ci_lwr,
      apc_ci_upr
    )

  aapc <- estimate_aapc(
    apc_data,
    deg_free = deg_free,
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
      deg_free = deg_free,
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

segmented_reg_fit <- function(y, x, k, data, opts) {
  fit <- stats::lm(
    formula = rlang::inject(opts$fun(!!y) ~ lspline::lspline(!!x, k)),
    data = data
  )
  structure(fit, class = c("edgy_spline_fit", class(fit)))
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

#' @rdname segmented_reg
#' @export
extract_fits <- function(x, ...) {
  rlang::check_dots_empty()
  x$fits
}

#' @rdname segmented_reg
#' @export
extract_best_model <- function(x, ...) {
  rlang::check_dots_empty()
  x$best_fit$model
}

#' @rdname segmented_reg
#' @export
extract_best_predictions <- function(x, ...) {
  rlang::check_dots_empty()

  y_var <- rlang::f_lhs(x$formula)
  x_var <- rlang::f_rhs(x$formula)

  extract_best_model(x) |>
    broom::augment() |>
    dplyr::transmute(
      !!x_var := x$x_vals,
      !!y_var := x$y_vals,
      .pred = x$opts$inv_fun(.fitted)
    )
}


#' @rdname segmented_reg
#' @export
extract_best_metrics <- function(x, ...) {
  rlang::check_dots_empty()
  x$best_fit$metrics
}

#' @rdname segmented_reg
#' @export
extract_best_apc <- function(x, ...) {
  rlang::check_dots_empty()
  x$best_fit$APC
}

#' @rdname segmented_reg
#' @export
extract_best_aapc <- function(x, ...) {
  rlang::check_dots_empty()
  x$best_fit$AAPC
}
