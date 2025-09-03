#' Segmented Regression using Linear Splines
#'
#' @param formula Model formula.
#' @param data Data frame or tibble.
#' @param opts A `knots_opt()` object providing Knot location options.
#' @param metric Metric to use for model selection. One of `"bic"`, `"bic3"`,
#'   `"wbic"`, `"aicc"`, or `"aicc".` Default is `"bic"`.
#' @param conf_level The confidence level to use for the confidence interval.
#'   Must be strictly greater than 0 and less than 1.
#'   Defaults to 0.95 which corresponds to a 95 percent confidence interval.
#' @param log_y Logical for whether to fit `log(y) ~ x` instead of `y ~ x`.
#'   Default is `TRUE`.
#' @param x A `segmented_reg()` object.
#' @param ... Not currently used.
#'
#' @return A named-list of class `"edgy_segmented_reg"`.
#'
#' @importFrom rlang :=
#'
#' @md
#' @export
segmented_reg <- function(
  formula,
  data,
  opts = knot_opts(),
  metric = 'bic',
  conf_level = 0.95,
  log_y = TRUE,
  ...
) {
  x_var <- rlang::f_rhs(formula)
  x_vals <- rlang::inject(`$`(data, !!x_var))
  y_var <- rlang::f_lhs(formula)
  y_vals <- rlang::inject(`$`(data, !!y_var))

  if (log_y) {
    opts$inv_fun <- exp
    opts$fun <- log
  } else {
    opts$inv_fun <- identity
    opts$fun <- identity
  }

  opts$conf.type <- "parametric"

  knot_outputs <- make_knot_sets(x_vals, opts)
  opts <- knot_outputs$opts
  knot_sets <- knot_outputs$knot_sets

  no_knot_form <- rlang::inject(opts$fun(!!y_var) ~ !!x_var)

  no_knot_model <-
    tibble::enframe(list(stats::lm(no_knot_form, data)), NULL, 'model') |>
    dplyr::mutate(nknots = 0)

  res <-
    knot_sets |>
    tibble::enframe(NULL, 'knots') |>
    dplyr::mutate(
      model = purrr::map(knots, \(k) {
        segmented_reg_fit(y_var, x_var, k, data, opts)
      }),
      nknots = purrr::map_int(knots, length)
    ) |>
    dplyr::bind_rows(no_knot_model) |>
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
      bic = (k * log(N)) - (2 * L),
      bic3 = log(sse / N) + ((3 * k + 2) / N) * log(N),
      w = 0.5,
      wbic = (bic * (1 - w)) + (bic3 * w)
    )

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
    periods = periods,
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
      apc_ci_upr,
      apc_pval
    )

  aapc <- estimate_aapc(
    apc_data,
    deg_free = deg_free,
    conf_level = conf_level,
    opts = opts
  )

  ret <- list(
    formula = formula,
    data = data_with_preds,
    y = y_vals,
    x = x_vals,
    opts = opts,
    knot_sets = knot_sets,
    fits = res,
    best_fit = list(
      fit = best_fit,
      model = best_fit$model[[1]],
      criterion = metric,
      conf_level = conf_level,
      deg_free = deg_free,
      metrics = list(
        nknots = best_fit$nknots,
        knots = best_fit$knots[[1]],
        nobs = best_fit$N,
        logLik = best_fit$L,
        AIC = best_fit$aic,
        AICc = best_fit$aicc,
        BIC = best_fit$bic,
        BIC3 = best_fit$bic3,
        WBIC = best_fit$wbic
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
    rlang::inject(opts$fun(!!y) ~ lspline::lspline(!!x, k)),
    data
  )
  structure(fit, class = c("edgy_spline_fit", class(fit)))
}

#' @rdname segmented_reg
#' @export
print.edgy_segmented_reg <- function(x, ...) {
  nmodels <- format(nrow(x$fits), big.mark = ',')
  mets <- x$best_fit$metrics
  locs <- glue::glue_collapse(mets$knots, ", ")
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
    broom::augment(se = TRUE) |>
    dplyr::transmute(
      !!y_var := x$y,
      !!x_var := x$x,
      .fitted = x$opts$inv_fun(.fitted),
      .fitted_lci = .fitted - stats::qnorm(0.975) * x$opts$inv_fun(.se.fit),
      .fitted_uci = .fitted + stats::qnorm(0.975) * x$opts$inv_fun(.se.fit)
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
