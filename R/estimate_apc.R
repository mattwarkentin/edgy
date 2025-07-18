estimate_apc <- function(x, knots, model, periods, df, conf.level) {
  periods <- get_periods(x, knots)

  broom::tidy(model) |>
    dplyr::filter(term != '(Intercept)') |>
    dplyr::bind_cols(periods) |>
    dplyr::mutate(
      var = diag(stats::vcov(model))[-1],
      apc = (exp(estimate) - 1) * 100,
      apc_ci_lwr = (exp(
        (estimate - (stats::qt(1 - (1 - conf.level) / 2, df) * std.error))
      ) -
        1) *
        100,
      apc_ci_upr = (exp(
        (estimate + (stats::qt(1 - (1 - conf.level) / 2, df) * std.error))
      ) -
        1) *
        100
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
