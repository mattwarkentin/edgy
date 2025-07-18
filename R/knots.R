#' Knot Options
#'
#' Define rules for selecting knot locations for segmented regression.
#'
#' @param min_knots Minimum number of knots to assess. Default is `1L`.
#' @param max_knots Maximum number of knots to assess. If `NULL`, the maximum
#'   number of knots will be estimated based on the number of data points.
#' @param min_obs_end Minimum number of points from either end of the observed
#'   data range to allow a knot location. Default is `2L`.
#' @param min_obs_between Number of observed data points between knot locations.
#'   Default is `2L`.
#' @param pts_between Number of points to place between adjacent observed values.
#'   Default is `0L`.
#' @param x Observed time points.
#' @param opts A `knot_opts()` object.
#' @param force Force the knot locations to be generated despite guidelines
#'   around the mininum number of observaions recommended for analysis.
#'   Default is `FALSE`.
#' @param ... Not currently used.
#'
#' @details
#' The minimum number of observations required to run the segmented regression
#'   analysis is estimated as `(2 * min_obs_end) + (max_knots - 1) * min_obs_between + max_knots`.
#'
#' @return `knot_opts()` returns a list with class `"edgy_knot_opts"`.
#'   `make_knot_sets()` returns a list of integer vectors that specify the
#'   knot locations to evaluate in grid search.
#'
#' @md
#'
#' @export
knot_opts <- function(
  min_knots = 1L,
  max_knots = NULL,
  min_obs_end = 2L,
  min_obs_between = 2L,
  pts_between = 0L
) {
  if (!rlang::is_null(max_knots)) {
    check_scalar_integer(max_knots, 'max_knots')
  }
  check_scalar_integer(min_knots, 'min_knots')
  check_scalar_integer(min_obs_end, 'min_obs_end')
  check_scalar_integer(min_obs_between, 'min_obs_between')
  check_scalar_integer(pts_between, 'pts_between')

  if (min_knots < 1L) {
    rlang::abort(glue::glue('`min_knots` must be greater than or equal to 1'))
  }

  if (!rlang::is_null(max_knots) && max_knots <= min_knots) {
    rlang::abort(glue::glue('`max_knots` must be greater than `min_knots`'))
  }

  if (min_obs_end < 0) {
    rlang::abort(glue::glue(
      '`min_obs_end` must be greater than or equal to zero'
    ))
  }

  if (min_obs_end < 0) {
    rlang::abort(glue::glue(
      '`min_obs_end` must be greater than or equal to zero'
    ))
  }

  if (pts_between < 0L) {
    rlang::abort(glue::glue('`pts_between` must be greater than or equal to 1'))
  }

  structure(
    list(
      min_knots = min_knots,
      max_knots = max_knots,
      min_obs_end = min_obs_end,
      min_obs_between = min_obs_between,
      pts_between = pts_between
    ),
    class = "edgy_knot_opts"
  )
}

#' @rdname knot_opts
#' @export
print.edgy_knot_opts <- function(x, ...) {
  if (rlang::is_null(x$max_knots)) {
    max_knots <- "<<Default>>"
  }

  cli::cli_rule('Segmented Regression Knot Options')
  cli::cat_bullet(glue::glue('Minimum number of knots: {x$min_knots}'))
  cli::cat_bullet(glue::glue('Maximum number of knots: {max_knots}'))
  cli::cat_bullet(glue::glue(
    'Minimum observations from end points: {x$min_obs_end}'
  ))
  cli::cat_bullet(glue::glue(
    'Minimum observations between knots: {x$min_obs_between}'
  ))
  cli::cat_bullet(glue::glue('Points between adjacent values: {x$pts_between}'))
}

#' @rdname knot_opts
#' @export
make_knot_sets <- function(x, opts, force = FALSE) {
  rlang::inherits_any(opts, "edgy_knot_opts")

  suggested_max_knots <- suggest_max_knots(x)

  if (rlang::is_null(opts$max_knots)) {
    opts$max_knots <- suggested_max_knots
  }

  if (opts$max_knots <= opts$min_knots) {
    rlang::abort(glue::glue('`max_knots` must be greater than `min_knots`'))
  }

  if (opts$max_knots > suggested_max_knots) {
    rlang::warn(glue::glue(
      "`max_knots` ({opts$max_knots}) is larger than the recommended maximum number of knots ({suggested_max_knots})"
    ))
  }

  min_num_obs <-
    (2 * opts$min_obs_end) +
    (opts$max_knots - 1) * opts$min_obs_between +
    opts$max_knots

  if (length(x) < min_num_obs & !force) {
    rlang::abort(
      'Ther are too few observations to run segmented regression with current knot options'
    )
  }

  knot_set <- seq(opts$min_knot, opts$max_knot, by = 1)

  knot_sets <-
    purrr::map(knot_set, \(k) combinat::combn(x, k, simplify = FALSE)) |>
    purrr::list_flatten() |>
    purrr::keep(\(k) {
      all(
        k >= min(x) + opts$min_obs_end &
          k <= max(x) - opts$min_obs_end
      )
    }) |>
    purrr::keep(\(k) {
      if (base::length(k) == 1) {
        return(TRUE)
      }
      lk <- base::length(k) - 1
      for (i in 1:lk) {
        if (k[[i + 1]] - k[[i]] < opts$min_obs_between) return(FALSE)
      }
      TRUE
    })
  knot_sets
}

suggest_max_knots <- function(x) {
  nobs <- length(x)
  dplyr::case_when(
    dplyr::between(nobs, 0, 6) ~ 0,
    dplyr::between(nobs, 7, 11) ~ 1,
    dplyr::between(nobs, 12, 16) ~ 2,
    dplyr::between(nobs, 17, 21) ~ 3,
    dplyr::between(nobs, 22, 26) ~ 4,
    dplyr::between(nobs, 27, 31) ~ 5,
    dplyr::between(nobs, 32, 36) ~ 6,
    dplyr::between(nobs, 37, Inf) ~ 7
  )
}
