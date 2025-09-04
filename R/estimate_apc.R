estimate_apc <- function(x, knots, model, deg_free, conf_level, opts) {
  periods <- get_periods(x, knots)

  broom::tidy(model) |>
    dplyr::filter(term != '(Intercept)') |>
    dplyr::bind_cols(periods) |>
    dplyr::mutate(
      var = diag(stats::vcov(model))[-1],
      apc = (opts$inv_fun(estimate) - 1) * 100,
      apc_se = opts$inv_fun(estimate) * std.error,
      apc_ci_lwr = (opts$inv_fun(
        (estimate - (stats::qt(1 - (1 - conf_level) / 2, deg_free) * apc_se))
      ) -
        1) *
        100,
      apc_ci_upr = (opts$inv_fun(
        (estimate + (stats::qt(1 - (1 - conf_level) / 2, deg_free) * apc_se))
      ) -
        1) *
        100,
      apc_pval = NA_real_ #stats::dt(estimate, deg_free)
    )
}

get_periods <- function(x, knots) {
  period_years <- c(min(x), knots, max(x))

  tibble::tibble(period_years) |>
    dplyr::transmute(
      period_start = period_years,
      period_end = dplyr::lead(period_start)
    ) |>
    tidyr::drop_na()
}
