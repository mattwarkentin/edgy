#' Extract Information from Segmented Regression Objects
#'
#' @param x An object with class `"egdy_segmented_reg"`
#' @param metric Which metric to use for selecting the optimal model fit?
#'   Default is `NULL` which uses the metric chosen in `segmented_reg()`.
#' @param nknots Number of knots used in model fit. Default is `NULL` which
#'   returns the globally optimal model. Set this to a specific value to extract
#'   the optimal model for a given number of knots.
#' @param ... Not currently used.
#'
#' @return The following objects are returned by the extraction functions:
#'   - `extract_fits()`: a `tibble` of statistics for model fits
#'   - `extract_best_fit()`: a one-row `tibble` of statistics for a specific
#'     model fit
#'   - `extract_best_model*()`: a model object with class `"cantrends_spline_fit"`
#'   - `extract_best_metrics()`: a named-list of model fit metrics
#'   - `extract_best_predictions()`: a `tibble` with the variables used in the
#'     model formula and predictions based on the fitted model
#'   - `extract_best_apc()`: a `tibble` with data for the Annual Percent Change
#'     (APC) for the chosen model
#'   - `extract_best_aapc()`: a `tibble` with data for the Avergae Annual
#'     Percent Change (AAPC) for the chosen model
#'
#' @export
#'
#' @examples
#' df <- read.delim(system.file("example.txt", package = "cantrends"), header = FALSE)
#' res <- segmented_reg(V3 ~ V2, data = df)
#'
#' extract_best_fit(res)
#' extract_best_model(res)
#' extract_best_metrics(res)
#' extract_best_predictions(res)
#' extract_best_apc(res)
#' extract_best_aapc(res)
extract_fits <- function(x, ...) {
  rlang::check_dots_empty()
  x$fits
}

#' @rdname extract_fits
#' @export
extract_best_fit <- function(x, metric = NULL, nknots = NULL, ...) {
  rlang::check_dots_empty()

  if (rlang::is_null(metric) & rlang::is_null(nknots)) {
    return(x$best_fit$fit)
  }

  if (rlang::is_null(metric)) {
    metric <- x$best_fit$criterion
  }

  if (rlang::is_null(nknots)) {
    nknots <- x$best_fit$metrics$nknots
  }

  if (nknots > x$opts$max_knots) {
    rlang::abort(
      "`nknots` is larger than `opts$max_knots` used during fitting."
    )
  }

  if (!metric %in% x$opts$metrics) {
    metrics <- glue::glue_collapse(x$opts$metrics, ", ")
    rlang::abort(
      message = glue::glue("`metric` must be one of {metrics}")
    )
  }

  x$fits |>
    dplyr::filter(.data$nknots == .env$nknots) |>
    dplyr::slice_min(!!rlang::sym(metric), with_ties = FALSE)
}

#' @rdname extract_fits
#' @export
extract_best_model <- function(x, metric = NULL, nknots = NULL, ...) {
  rlang::check_dots_empty()

  if (rlang::is_null(metric) & rlang::is_null(nknots)) {
    return(x$best_fit$model)
  }

  extract_best_fit(x, metric, nknots)$model[[1]]
}

#' @rdname extract_fits
#' @export
extract_best_metrics <- function(x, metric = NULL, nknots = NULL, ...) {
  rlang::check_dots_empty()

  if (rlang::is_null(metric) & rlang::is_null(nknots)) {
    return(x$best_fit$metrics)
  }

  curr_fit <- extract_best_fit(x, metric, nknots)

  list(
    nknots = curr_fit$nknots,
    nobs = curr_fit$N,
    logLik = curr_fit$L,
    AIC = curr_fit$aic,
    AICc = curr_fit$aicc,
    BIC = curr_fit$bic,
    BIC3 = curr_fit$bic3
  )
}

#' @rdname extract_fits
#' @export
extract_best_predictions <- function(x, metric = NULL, nknots = NULL, ...) {
  rlang::check_dots_empty()

  y_var <- rlang::f_lhs(x$formula)
  x_var <- rlang::f_rhs(x$formula)

  curr_fit <- extract_best_fit(x, metric, nknots)
  periods <- get_periods(x$x_vals, curr_fit$knots[[1]])

  extract_best_model(x, metric, nknots) |>
    broom::augment() |>
    dplyr::transmute(
      !!x_var := x$x_vals,
      !!y_var := x$y_vals,
      .pred = x$opts$inv_fun(.fitted)
    ) |>
    dplyr::left_join(
      periods,
      by = dplyr::join_by(dplyr::between(!!x_var, period_start, period_end))
    )
}

#' @rdname extract_fits
#' @export
extract_best_apc <- function(x, metric = NULL, nknots = NULL, ...) {
  rlang::check_dots_empty()

  if (rlang::is_null(metric) & rlang::is_null(nknots)) {
    return(x$best_fit$APC)
  }

  curr_fit <- extract_best_fit(x, metric, nknots)

  estimate_apc(
    x = x$x_vals,
    knots = curr_fit$knots[[1]],
    model = curr_fit$model[[1]],
    deg_free = curr_fit$df,
    conf_level = x$best_fit$conf_level,
    opts = x$opts
  ) |>
    dplyr::select(
      segment,
      period_start,
      period_end,
      apc,
      apc_ci_lwr,
      apc_ci_upr,
      apc_pval
    )
}

#' @rdname extract_fits
#' @export
extract_best_aapc <- function(x, metric = NULL, nknots = NULL, ...) {
  rlang::check_dots_empty()

  if (rlang::is_null(metric) & rlang::is_null(nknots)) {
    return(x$best_fit$AAPC)
  }

  curr_fit <- extract_best_fit(x, metric, nknots)

  apc_data <- estimate_apc(
    x = x$x_vals,
    knots = curr_fit$knots[[1]],
    model = curr_fit$model[[1]],
    deg_free = curr_fit$df,
    conf_level = x$best_fit$conf_level,
    opts = x$opts
  )

  estimate_aapc(
    apc_data,
    deg_free = curr_fit$df,
    conf_level = x$best_fit$conf_level,
    opts = x$opts
  )
}
