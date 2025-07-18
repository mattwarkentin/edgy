#' Define knot placements
#'
#' @param min_knot Fewest number of knots to assess. Default is `0L`.
#' @param max_knot Largest number of knots to assess. Default is `4L`.
#' @param min_obs_end Minimum number of points from either end of the observed
#'   data range to allow a knot location. Default is `2L`.
#' @param min_obs_between Number of observed data points between knot locations.
#'   Default is `2L`.
#' @param pts_between Number of theoretical points between knot locations.
#'   Default is `4L`.
#' @param x Time points.
#' @param opts A `knot_opts()` object.
#'
#' @md
#'
#' @export
knot_opts <- function(
  min_knot = 0L,
  max_knot = 4L,
  min_obs_end = 2L,
  min_obs_between = 2L,
  pts_between = 0L
) {
  structure(
    list(
      min_knot = check_scalar_integer(min_knot, 'min_knot'),
      max_knot = check_scalar_integer(max_knot, 'max_knot'),
      min_obs_end = check_scalar_integer(min_obs_end, 'min_obs_end'),
      min_obs_between = check_scalar_integer(
        min_obs_between,
        'min_obs_between'
      ),
      pts_between = check_scalar_integer(pts_between, 'pts_between')
    ),
    class = "edgy_knot_opts"
  )
}

#' @rdname knot_opts
#' @export
make_knot_grid <- function(x, opts) {
  rlang::inherits_any(opts, "edgy_knot_opts")

  knot_set <- seq(opts$min_knot, opts$max_knot, by = 1)

  knot_grid <-
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
  knot_grid
}
