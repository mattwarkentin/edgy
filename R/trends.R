#' Trend Regression using Linear Splines
#'
#' @param formula Model formula.
#' @param data Data frame or tibble.
#' @param opts A `knots_opt()` object providing Knot location options.
#' @param metric Metric to use for model selection. Default is `"bic"`. One
#'   of `"bic"`, `"aicc"`, or `"aicc".`
#' @param ... Not currently used.
#'
#' @return A named-list with components `metrics`, `data`, and `APC`.
#'
#' @md
#' @export
trend_reg <- function(
    formula,
    data,
    metric = 'bic',
    opts = knot_opts(),
    ...
) {
  x_var <- rlang::f_rhs(formula)
  x_vals <- rlang::inject(`$`(data, !!x_var))
  y_var <- rlang::f_lhs(formula)

  opts$inv_fun <- identity
  if (stringr::str_detect(rlang::expr_text(y_var), '^log')) opts$inv_fun <- exp

  knot_set <- seq(opts$min_knot, opts$max_knot, by = 1)

  knot_grid <-
    purrr::map(knot_set, \(k) combinat::combn(x_vals, k, simplify = FALSE)) |>
    purrr::list_flatten() |>
    purrr::keep(\(k) {
      all(k >= min(x_vals) + opts$min_obs_end &
            k <= max(x_vals) - opts$min_obs_end)
    }) |>
    purrr::keep(\(k) {
      if (base::length(k) == 1) return(TRUE)
      lk <- base::length(k) - 1
      for (i in 1:lk) {
        if (k[[i + 1]] - k[[i]] < opts$min_obs_between) return(FALSE)
      }
      TRUE
    })

  no_knot_form <- rlang::inject(!!y_var ~ !!x_var)

  no_knot_model <-
    tibble::enframe(list(stats::lm(no_knot_form, data)), NULL, 'model') |>
    dplyr::mutate(nknots = 0)

  res <-
    knot_grid |>
    tibble::enframe(NULL, 'knots') |>
    dplyr::mutate(
      model = purrr::map(knots, \(k) {
        stats::lm(rlang::inject(!!y_var ~ lspline::lspline(!!x_var, k)), data)
      }),
      nknots = purrr::map_int(knots, length)
    ) |>
    dplyr::bind_rows(no_knot_model) |>
    dplyr::mutate(
      preds = purrr::map(model, \(m) suppressWarnings(stats::predict(m, interval = 'prediction'))),
      preds = purrr::map(preds, tibble::as_tibble),
      k = purrr::map_int(model, \(m) base::length(m$coefficients) + 1),
      L = purrr::map_dbl(model, \(m) as.numeric(stats::logLik(m))),
      N = purrr::map_int(model, stats::nobs),
      aic = -2 * L + 2 * k,
      aicc = aic + 2 * (k * (k + 1)) / (N - k - 1),
      bic = (k * log(N)) - (2 * L)
    )

  best_fit <- dplyr::slice_min(res, bic)

  data_with_preds <-
    dplyr::bind_cols(
      data,
      best_fit |>
        dplyr::select(preds) |>
        tidyr::unnest(preds)
    ) |>
    dplyr::mutate(dplyr::across(c(fit, lwr, upr), opts$inv_fun))

  period_years <- c(min(x_vals), best_fit$knots[[1]], max(x_vals))

  periods <-
    tibble::tibble(period_years) |>
    dplyr::transmute(
      period_start = period_years,
      period_end = dplyr::lead(period_start)
    ) |>
    tidyr::drop_na()

  df <- df.residual(best_fit$model[[1]])

  apc_data <-
    broom::tidy(best_fit$model[[1]]) |>
    dplyr::filter(term != '(Intercept)') |>
    dplyr::bind_cols(periods) |>
    dplyr::mutate(
      apc = (exp(estimate) - 1) * 100,
      apc_lwr = (exp((estimate - (qt(0.975, df) * std.error))) -1) * 100,
      apc_upr = (exp((estimate + (qt(0.975, df) * std.error))) -1) * 100
    )

  apc <-
    apc_data |>
    dplyr::select(period_start, period_end, apc, apc_lwr, apc_upr)

  aapc <-
    apc_data |>
    dplyr::mutate(
      d = period_end - period_start,
      db = d * estimate
    ) |>
    dplyr::summarise(
      period_start = min(period_start),
      period_end = max(period_end),
      aapc = (exp(sum(db) / (sum(d))) - 1) * 100
    )

  list(
    metrics = list(
      knots = best_fit$knots[[1]],
      nobs = best_fit$N,
      logLik = best_fit$L,
      AIC = best_fit$aic,
      AICc = best_fit$aicc,
      BIC = best_fit$bic
    ),
    model = best_fit$model[[1]],
    data = data_with_preds,
    APC = apc,
    AAPC = aapc
  )
}

globalVariables(c('min_knots', 'max_knots', 'knots', 'model', 'preds',
                  'L', 'N', 'aic', 'bic', 'apc', 'change', 'fit',
                  'k', 'left_point', 'right_point', 'lwr', 'upr',
                  'nyears', 'period_end'))

#' Spline knot options
#'
#' @param min_knot Fewest number of knots to assess. Default is `0L`.
#' @param max_knot Largest number of knots to assess. Default is `4L`.
#' @param min_obs_end Minimum number of points from either end of the observed
#'   data range to allow a knot location. Default is `2L`.
#' @param min_obs_between Number of observed data points between knot locations.
#'   Default is `2L`.
#' @param pts_between Number of theoretical points between knot locations.
#'   Default is `4L`.
#'
#' @export
knot_opts <- function(
    min_knot = 0L,
    max_knot = 4L,
    min_obs_end = 2L,
    min_obs_between = 2L,
    pts_between = 0L
) {
  structure(
    list(
      min_knot = check_scalar_integer(min_knot, 'min_knot'),
      max_knot = check_scalar_integer(max_knot, 'max_knot'),
      min_obs_end = check_scalar_integer(min_obs_end, 'min_obs_end'),
      min_obs_between = check_scalar_integer(min_obs_between, 'min_obs_between'),
      pts_between = check_scalar_integer(pts_between, 'pts_between')
    ),
    class = 'edgy_knot_opts'
  )
}

check_scalar_integer <- function(x, name) {
  if (rlang::is_scalar_integerish(x)) return(x)
  rlang::abort(glue::glue('`{name}` must be a scalar integer.'))
}
