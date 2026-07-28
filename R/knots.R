#' Knot Options
#'
#' Define rules for selecting knot locations for segmented regression.
#'
#' @param min_knots Minimum number of knots to evaluate. Default is `1L`.
#' @param max_knots Maximum number of knots to evaluate. If `NULL`, the maximum
#'   number of knots will be estimated based on the number of data points.
#' @param min_obs_end Minimum number of points from either end of the observed
#'   data range to allow a knot location. Default is `2L`.
#' @param min_obs_between Minimum number of data points between knot
#'   locations. Default is `2L`.
#' @param pts_between Number of new data points to place between adjacent
#'   observed data points. Default is `0L`.
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
#'   `make_knot_sets()` returns a two-item list with (1) a list of vectors that
#'   specify the knot locations to evaluate in grid search and (2) the updated
#'   `opts` object. `make_point_set()` returns the vector of points used to
#'   generate the knot sets`.
#'
#' @md
#'
#' @examples
#' knot_opts()
#'
#' @export
knot_opts <- function(
  min_knots = 1L,
  max_knots = NULL,
  min_obs_end = 2L,
  min_obs_between = 2L,
  pts_between = 0L,
  force = FALSE
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

  if (!rlang::is_null(max_knots) && max_knots < min_knots) {
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
    rlang::abort(glue::glue('`pts_between` must be greater than or equal to 0'))
  }

  structure(
    list(
      min_knots = min_knots,
      max_knots = max_knots,
      min_obs_end = min_obs_end,
      min_obs_between = min_obs_between,
      pts_between = pts_between,
      force = force
    ),
    class = "edgy_knot_opts"
  )
}

#' @rdname knot_opts
#' @export
print.edgy_knot_opts <- function(x, ...) {
  if (rlang::is_null(x$max_knots)) {
    max_knots <- "<<Default>>"
  } else {
    max_knots <- x$max_knots
  }

  cli::cli_rule("Segmented Regression Knot Options")
  cli::cat_bullet(glue::glue("Minimum number of knots: {x$min_knots}"))
  cli::cat_bullet(glue::glue("Maximum number of knots: {max_knots}"))
  cli::cat_bullet(glue::glue(
    "Minimum observations from end points: {x$min_obs_end}"
  ))
  cli::cat_bullet(glue::glue(
    "Minimum observations between knots: {x$min_obs_between}"
  ))
  cli::cat_bullet(glue::glue("Points between adjacent values: {x$pts_between}"))
}

#' @rdname knot_opts
#' @export
make_knot_sets <- function(x, opts) {
  rlang::inherits_any(opts, "edgy_knot_opts")

  suggested_max_knots <- suggest_max_knots(x)

  if (rlang::is_null(opts$max_knots)) {
    opts$max_knots <- suggested_max_knots
  }

  if (opts$max_knots < opts$min_knots) {
    rlang::abort(glue::glue("`max_knots` must be greater than `min_knots`."))
  }

  if (opts$max_knots > suggested_max_knots) {
    rlang::warn(glue::glue(
      "`max_knots` ({opts$max_knots}) is larger than the recommended maximum number of knots ({suggested_max_knots})."
    ))
  }

  min_num_obs <-
    (2 * opts$min_obs_end) +
    (opts$max_knots - 1) * opts$min_obs_between +
    opts$max_knots

  if (length(x) < min_num_obs & !opts$force) {
    rlang::abort(
      "There are too few observations to run segmented regression with current knot options. Set `knot_opts(force = TRUE)` to override this behavior."
    )
  }

  num_knots <- seq(opts$min_knots, opts$max_knots, by = 1)

  point_set <- make_point_set(x, opts)

  knot_sets <- purrr::map(
    .x = num_knots,
    .f = \(k) make_knots_for_k(point_set, k, opts)
  ) |>
    purrr::list_flatten()

  list(knot_sets = knot_sets, opts = opts)
}

make_knots_for_k <- function(x, k, opts) {
  n <- length(x)

  # interior length
  m <- n - 2L * opts$min_obs_end

  if (m <= 0L || k <= 0L) {
    return(list())
  }

  # compressed length
  max_base <- m - (opts$min_obs_between - 1L) * (k - 1L)

  if (max_base < k) {
    return(list())
  }

  base_combos <- utils::combn(seq_len(max_base), k, simplify = FALSE)

  result <- purrr::map(
    .x = base_combos,
    .f = \(y) {
      expanded <- y + (opts$min_obs_between - 1L) * (seq_along(y) - 1L)
      idx <- as.integer(expanded + opts$min_obs_end)
      x[idx]
    }
  )

  result
}

#' @rdname knot_opts
#' @export
make_point_set <- function(x, opts) {
  unique_obs_pts <- sort(unique(x))

  point_set <- unique_obs_pts

  if (opts$pts_between > 0) {
    seq_pairs <- purrr::map2(utils::head(x, -1), utils::tail(x, -1), c)

    new_pts <-
      purrr::map(seq_pairs, \(x) {
        left <- x[[1]]
        right <- x[[2]]
        spacing <- right - left
        seq(left, right, length.out = 2 + opts$pts_between)
      })

    point_set <- purrr::list_c(new_pts)
  }

  point_set
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
