estimate_aapc <- function(data, deg_free, conf_level, opts) {
  z <- 1 - (1 - conf_level) / 2

  # In the Joinpoint software, the AAPC confidence interval is based on the
  # normal distribution, and the APC confidence interval is based on a
  # t distribution. If an AAPC lies entirely within a single joinpoint segment,
  # the AAPC is equal to the APC for that segment. To obtain consistency between
  # the APC and AAPC confidence intervals in this situation, the confidence
  # interval for the AAPC has been modified to be identical to that used for
  # the APC using the t distribution instead of the normal distribution.
  num_segments <- nrow(data)

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
      aapc_ci_lwr = ifelse(
        num_segments == 1L,
        apc_ci_lwr,
        (opts$inv_fun(
          opts$fun((aapc / 100) + 1) - (stats::qnorm(z) * aapc_se)
        ) -
          1) *
          100
      ),
      aapc_ci_upr = ifelse(
        num_segments == 1L,
        apc_ci_upr,
        (opts$inv_fun(
          opts$fun((aapc / 100) + 1) + (stats::qnorm(z) * aapc_se)
        ) -
          1) *
          100
      )
    ) |>
    dplyr::select(-aapc_se)
}
