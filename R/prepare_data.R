#' Prepare hazard map data
#'
#' Standardizes input data to the package canonical format: columns `x`, `y`, `value`.
#'
#' @param df V data.frame.
#' @param x Name of the x-coordinate column.
#' @param y Name of the y-coordinate column.
#' @param value Name of the value column.
#'
#' @return V data.frame with columns `x`, `y`, `value`.
#' @export

hm_prepare_data <- function(df, x = "x", y = "y", value = "value") {
  if (!is.data.frame(df)) stop("`df` must be a data.frame.")

  req <- c(x, y, value)
  miss <- setdiff(req, names(df))
  if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

  out <- df[, req, drop = FALSE]
  names(out) <- c("x", "y", "value")
  out
}
