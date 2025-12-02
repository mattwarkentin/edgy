estimate_aapc <- function(data, deg_free, conf_level, opts) {
  z_val <- 1 - (1 - conf_level) / 2

  data |>
    dplyr::mutate(
      w = period_end - period_start,
      wn = w / sum(w),
      wb = w * estimate
    ) |>
    dplyr::summarise(
      period_start = min(period_start),
      period_end = max(period_end),
      aapc = (opts$inv_fun(sum(wb) / (sum(w))) - 1) * 100,
      aapc_se = sqrt(sum(wn^2 * beta_var)),
      aapc_ci_lwr = (opts$inv_fun(
        opts$fun((aapc / 100) + 1) - (stats::qnorm(z_val) * aapc_se)
      ) -
        1) *
        100,
      aapc_ci_upr = (opts$inv_fun(
        opts$fun((aapc / 100) + 1) + (stats::qnorm(z_val) * aapc_se)
      ) -
        1) *
        100
    )
}
