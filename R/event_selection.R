.hm_check_hazard <- function(x) {
  if (is.null(x) || !"hm_hazard" %in% class(x))
    stop("`x` must be an object of class `hm_hazard`.")
  if (is.null(x$data) || !is.array(x$data))
    stop("`x$data` must be an array.")
  if (is.null(x$coords) || !is.data.frame(x$coords))
    stop("`x$coords` must be a data.frame. Run `hm_standardize_coords()` first.")
  if (!inherits(x$time, c("Date", "POSIXct", "POSIXt")))
    stop("`x$time` must be decoded. Run `hm_decode_time()` first.")
  invisible(NULL)
}


.hm_to_flat_matrix <- function(x) {
  d     <- dim(x$data)
  T_len <- length(x$time)
  tidx  <- which(d == T_len)

  if (length(tidx) != 1L)
    stop("Cannot uniquely identify the time dimension. ",
         "Expected one dimension of length ", T_len,
         ", found: ", paste(d, collapse = " x "))

  sidx <- setdiff(seq_along(d), tidx)

  list(
    X    = matrix(as.vector(aperm(x$data, c(tidx, sidx))),
                  nrow = T_len,
                  ncol = prod(d[sidx])),
    tidx = tidx,
    sidx = sidx,
    d    = d
  )
}



#' Select extreme-event days by threshold
#'
#' Retains only the time steps in which at least one grid cell reaches or
#' exceeds \code{threshold}. Call after [hm_decode_time()] and before
#' [hm_to_matrix()].
#'
#' @param x An \code{hm_hazard} object after [hm_standardize_coords()] and
#'   [hm_decode_time()].
#' @param threshold A single finite numeric value.
#'
#' @return The same object with \code{$data} and \code{$time} restricted to
#'   extreme-event days. \code{$meta$filter} records the threshold and day
#'   counts for reproducibility.
#'
#' @seealso [hm_summarise_events()]
#' @export
hm_select_extreme_events <- function(x, threshold) {

  .hm_check_hazard(x)

  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold))
    stop("`threshold` must be a single finite numeric value.")

  flat <- .hm_to_flat_matrix(x)
  keep <- which(rowSums(flat$X >= threshold, na.rm = TRUE) > 0L)

  n_original <- length(x$time)
  n_retained <- length(keep)

  if (n_retained == 0L)
    stop("No days reach the threshold (", threshold, "). ",
         "Consider a lower value.")

  if (n_retained == n_original)
    message("All ", n_original, " days meet the threshold. ",
            "Filter recorded in $meta$filter.")

  idx                <- vector("list", length(flat$d))
  idx[flat$sidx]     <- lapply(flat$d[flat$sidx], seq_len)
  idx[[flat$tidx]]   <- keep

  x$data        <- do.call(`[`, c(list(x$data), idx, list(drop = FALSE)))
  x$time        <- x$time[keep]
  x$meta$filter <- list(
    threshold       = threshold,
    n_days_original = n_original,
    n_days_retained = n_retained
  )

  x
}



#' Summarise retained extreme-event days
#'
#' Returns a data.frame with one row per time step: date, number of cells
#' that met the threshold, and domain-wide maximum value.
#'
#' @param x An \code{hm_hazard} object, typically from
#'   [hm_select_extreme_events()].
#' @param threshold Threshold used to count exceeding cells. If \code{NULL},
#'   read from \code{x$meta$filter$threshold}.
#'
#' @return A data.frame with columns \code{date}, \code{n_cells_exceeded},
#'   and \code{max_value}.
#'
#' @seealso [hm_select_extreme_events()]
#' @export
hm_summarise_events <- function(x, threshold = NULL) {

  .hm_check_hazard(x)

  if (is.null(threshold)) {
    threshold <- x$meta$filter$threshold
    if (is.null(threshold))
      stop("`threshold` not found in `x$meta$filter`. ",
           "Run `hm_select_extreme_events()` first, or pass it directly.")
  }

  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold))
    stop("`threshold` must be a single finite numeric value.")

  X                <- .hm_to_flat_matrix(x)$X
  n_cells_exceeded <- as.integer(rowSums(X >= threshold, na.rm = TRUE))
  max_value        <- apply(X, 1L, max, na.rm = TRUE)
  max_value[!is.finite(max_value)] <- NA_real_

  data.frame(
    date             = x$time,
    n_cells_exceeded = n_cells_exceeded,
    max_value        = max_value,
    stringsAsFactors = FALSE
  )
}
