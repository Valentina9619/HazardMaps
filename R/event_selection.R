# Internal helper: validate that x is a ready-to-use hm_hazard object,
# meaning coords are standardized and time is decoded.
.hm_check_hazard <- function(x, call = "this function") {
  if (is.null(x) || !"hm_hazard" %in% class(x))
    stop("`x` must be an object of class `hm_hazard`.")
  if (is.null(x$data) || !is.array(x$data))
    stop("`x$data` must be an array.")
  if (is.null(x$coords) || !is.data.frame(x$coords))
    stop("`x$coords` must be a data.frame. Run `hm_standardize_coords()` first.")
  if (is.null(x$time) ||
      !(inherits(x$time, "Date") ||
        inherits(x$time, "POSIXct") ||
        inherits(x$time, "POSIXt")))
    stop("`x$time` must be a decoded Date/POSIXct. Run `hm_decode_time()` first.")
  invisible(NULL)
}


# Internal helper: find the time dimension index and return the data as a
# (time x cells) matrix. Used by both functions in this file.
.hm_to_flat_matrix <- function(x) {
  dat   <- x$data
  d     <- dim(dat)
  T_len <- length(x$time)

  tidx <- which(d == T_len)
  if (length(tidx) != 1L)
    stop(
      "Cannot uniquely identify the time dimension in `x$data`. ",
      "Expected one dimension of length ", T_len,
      ", found: ", paste(d, collapse = " x ")
    )

  sidx          <- setdiff(seq_along(d), tidx)
  dat_reordered <- aperm(dat, c(tidx, sidx))

  list(
    X    = matrix(as.vector(dat_reordered), nrow = T_len, ncol = prod(d[sidx])),
    tidx = tidx,
    sidx = sidx,
    d    = d
  )
}



#' Select extreme-event days by threshold
#'
#' Subsets an `hm_hazard` object, retaining only the time steps (days) on
#' which at least one grid cell in the full spatial domain has a hazard value
#' greater than or equal to `threshold`.
#'
#' @details
#' The filter is applied over the **entire spatial domain**. A day is retained
#' when at least one grid cell reaches or exceeds `threshold`, ignoring `NA`
#' values. Days where every cell is either below the threshold or `NA` are
#' dropped.
#'
#' The returned object is a valid `hm_hazard` with the same structure as the
#' input. `$data` is sliced along its time dimension, `$time` is replaced by
#' the filtered dates, and `$meta$filter` records the threshold and day counts
#' for reproducibility.
#'
#' Call this function **after** [hm_standardize_coords()] and
#' [hm_decode_time()], and **before** [hm_to_matrix()].
#'
#' @param x An `hm_hazard` object after [hm_standardize_coords()] and
#'   [hm_decode_time()].
#' @param threshold A single finite numeric value. Days where every grid cell
#'   is strictly below this value are removed.
#'
#' @return The same `hm_hazard` object with `$data` and `$time` restricted to
#'   extreme-event days. `$meta$filter` is added with three fields:
#'   `threshold`, `n_days_original`, and `n_days_retained`.
#'
#' @seealso [hm_summarise_events()] to inspect the retained days as a table.
#'
#' @examples
#' \dontrun{
#' f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
#'                  package = "HazardMaps")
#' x <- hm_read_netcdf(f, var = "i10fg")
#' x <- hm_standardize_coords(x)
#' x <- hm_decode_time(x)
#'
#' x_extreme <- hm_select_extreme_events(x, threshold = 25)
#' x_extreme$meta$filter
#' }
#'
#' @export
hm_select_extreme_events <- function(x, threshold) {

  .hm_check_hazard(x)

  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold))
    stop("`threshold` must be a single finite numeric value.")

  # Flatten to (time x cells) to evaluate the threshold row-wise
  flat <- .hm_to_flat_matrix(x)
  X    <- flat$X

  # A day is kept when at least one cell reaches the threshold
  keep <- which(rowSums(X >= threshold, na.rm = TRUE) > 0L)

  n_original <- length(x$time)
  n_retained <- length(keep)

  if (n_retained == 0L)
    stop(
      "No days reach the threshold (", threshold, "). ",
      "Consider using a lower value."
    )

  if (n_retained == n_original)
    message(
      "All ", n_original, " days meet the threshold. ",
      "The object is returned unchanged (filter is still recorded in $meta$filter)."
    )

  # Slice $data along the time dimension, keeping spatial dims intact
  idx         <- vector("list", length(flat$d))
  idx[flat$sidx]   <- lapply(flat$d[flat$sidx], seq_len)
  idx[[flat$tidx]] <- keep

  x$data <- do.call(`[`, c(list(x$data), idx, list(drop = FALSE)))
  x$time <- x$time[keep]

  x$meta$filter <- list(
    threshold       = threshold,
    n_days_original = n_original,
    n_days_retained = n_retained
  )

  x
}



#' Summarise retained extreme-event days
#'
#' Returns a data.frame with one row per time step, reporting the date, the
#' number of grid cells that reached or exceeded the threshold, and the
#' domain-wide maximum hazard value for that day.
#'
#' @details
#' Intended to be called on an object already filtered by
#' [hm_select_extreme_events()], in which case `threshold` is read
#' automatically from `x$meta$filter$threshold`. It can also be called on an
#' unfiltered object by supplying `threshold` explicitly, which summarises
#' every time step in the dataset.
#'
#' @param x An `hm_hazard` object after [hm_standardize_coords()] and
#'   [hm_decode_time()]. Typically the output of [hm_select_extreme_events()].
#' @param threshold A single finite numeric value used to count exceeding
#'   cells. If `NULL` (the default), the value stored in
#'   `x$meta$filter$threshold` is used. An error is raised if neither is
#'   available.
#'
#' @return A data.frame with one row per time step and three columns:
#'   \describe{
#'     \item{`date`}{Date of the event (`Date` or `POSIXct`).}
#'     \item{`n_cells_exceeded`}{Number of grid cells with value >= threshold.}
#'     \item{`max_value`}{Domain-wide maximum hazard value on that day.}
#'   }
#'
#' @seealso [hm_select_extreme_events()]
#'
#' @examples
#' \dontrun{
#' f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
#'                  package = "HazardMaps")
#' x <- hm_read_netcdf(f, var = "i10fg")
#' x <- hm_standardize_coords(x)
#' x <- hm_decode_time(x)
#'
#' x_extreme <- hm_select_extreme_events(x, threshold = 25)
#' hm_summarise_events(x_extreme)
#' }
#'
#' @export
hm_summarise_events <- function(x, threshold = NULL) {

  .hm_check_hazard(x)

  # Resolve threshold: from argument or from filter metadata
  if (is.null(threshold)) {
    threshold <- x$meta$filter$threshold
    if (is.null(threshold))
      stop(
        "`threshold` was not supplied and is not recorded in `x$meta$filter`. ",
        "Run `hm_select_extreme_events()` first, or pass `threshold` directly."
      )
  }

  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold))
    stop("`threshold` must be a single finite numeric value.")

  # Flatten to (time x cells)
  X <- .hm_to_flat_matrix(x)$X

  n_cells_exceeded <- as.integer(rowSums(X >= threshold, na.rm = TRUE))
  max_value        <- apply(X, 1L, max, na.rm = TRUE)

  # Rows where every cell is NA produce -Inf from max(); replace with NA
  max_value[!is.finite(max_value)] <- NA_real_

  data.frame(
    date             = x$time,
    n_cells_exceeded = n_cells_exceeded,
    max_value        = max_value,
    stringsAsFactors = FALSE
  )
}
