estimate_aapc <- function(data, df, conf.level) {
  data |>
    dplyr::mutate(
      d = period_end - period_start,
      db = d * estimate
    ) |>
    dplyr::summarise(
      period_start = min(period_start),
      period_end = max(period_end),
      aapc = (exp(sum(db) / (sum(d))) - 1) * 100,
      aapc_ci_lwr = (exp(
        log((aapc / 100) + 1) -
          (stats::qt(1 - (1 - conf.level) / 2, df) * sqrt(sum(d^2 * var)))
      ) -
        1) *
        100,
      aapc_ci_upr = (exp(
        log((aapc / 100) + 1) +
          (stats::qt(1 - (1 - conf.level) / 2, df) * sqrt(sum(d^2 * var)))
      ) -
        1) *
        100
    )
}
