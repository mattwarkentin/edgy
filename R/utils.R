check_scalar_integer <- function(x, name) {
  if (rlang::is_scalar_integerish(x)) {
    return(x)
  }
  rlang::abort(glue::glue('`{name}` must be a scalar integer.'))
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
  'value',
  'name',
  '.fitted',
  'nknots'
))
