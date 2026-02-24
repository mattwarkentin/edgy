check_scalar_integer <- function(x, name) {
  if (rlang::is_scalar_integerish(x)) {
    return(x)
  }
  rlang::abort(glue::glue('`{name}` must be a scalar integer.'))
}

zero_rate_transform <- function(n, num = 0.5) {
  0.5 / n
}

check_positive_rates <- function(x) {
  if (rlang::is_double(x) & all(x > 0)) {
    return(TRUE)
  }
  rlang::abort(
    message = "Rates must be strictly positive. Consider transforming zero rates (e.g., 0.5 / n, where n is the number of persons in that year).",
    call = rlang::env_parent()
  )
}

globalVariables(c(
  'min_knots',
  'max_knots',
  'knots',
  'model',
  'preds',
  'L',
  'N',
  'aic',
  'bic',
  'apc',
  'change',
  'fit',
  'k',
  'left_point',
  'right_point',
  'lwr',
  'upr',
  'nyears',
  'period_end',
  'period_start',
  'd',
  'estimate',
  'db',
  'term',
  'std.error',
  'sse',
  'w',
  'bic3',
  'periods',
  'apc_ci_lwr',
  'apc_ci_upr',
  'apc_se',
  'apc_pval',
  'var',
  'aapc',
  'aapc_se',
  'beta_var',
  'wb',
  'wn',
  'value',
  'name',
  '.fitted',
  'nknots'
))
