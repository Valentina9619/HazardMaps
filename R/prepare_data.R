#' Read a hazard NetCDF dataset
#'
#' Reads a NetCDF file containing a hazard variable on a spatial grid
#' (regular or rotated) and returns a basic `hm_hazard` object.
#' CRS is NOT required.
#'
#' @details
#' The function supports NetCDF files following common CF conventions for
#' spatial coordinates. Latitude and longitude may be provided either as
#' coordinate variables or as dimension coordinates, and may be encoded as:
#' \itemize{
#'   \item one-dimensional regular grids (e.g., Copernicus / ERA5 products),
#'   \item two-dimensional curvilinear or rotated grids (e.g., EURO-CORDEX).
#' }
#' Coordinate names \code{lat}/\code{lon} and \code{latitude}/\code{longitude}
#' are automatically detected.
#'
#' @param file Path to a NetCDF file.
#' @param var Optional. Name of the hazard variable. If `NULL` and the file
#'   contains exactly one variable, it is selected automatically. If multiple
#'   variables exist and `var` is `NULL`, an error is raised.
#'
#' @return An object of class `hm_hazard`, a named list with elements:
#'   \describe{
#'     \item{`data`}{Raw array read from the NetCDF file.}
#'     \item{`coords`}{Named list with `lat`, `lon`, `rlat`, `rlon` (raw).}
#'     \item{`time`}{Named list with `values` and `units` (raw CF encoding).}
#'     \item{`meta`}{Named list with `variable`, `grid_type`, `reader`, `file`.}
#'   }
#'
#' @seealso [hm_standardize_coords()], [hm_decode_time()]
#' @export
hm_read_netcdf <- function(file, var = NULL) {

  if (!is.character(file) || length(file) != 1L)
    stop("`file` must be a single path string.")
  if (!file.exists(file))
    stop("File not found: ", file)

  nc <- ncdf4::nc_open(file)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  # --- variable selection -----------------------------------------------------
  vnames <- names(nc$var)

  if (is.null(var)) {
    if (length(vnames) != 1L)
      stop(
        "Multiple variables found. Please specify `var`. ",
        "Candidates: ", paste(vnames, collapse = ", ")
      )
    var <- vnames[[1L]]
  } else {
    if (!var %in% vnames)
      stop(
        "Variable `", var, "` not found. ",
        "Candidates: ", paste(vnames, collapse = ", ")
      )
  }

  # --- grid type detection ----------------------------------------------------
  dnames <- names(nc$dim)

  if (!"time" %in% dnames)
    stop("No `time` dimension found in the NetCDF file.")

  grid_type <- if (all(c("rlat", "rlon") %in% dnames)) "rotated" else "regular"

  # --- read time axis ---------------------------------------------------------
  time_vals  <- ncdf4::ncvar_get(nc, "time")
  time_units <- ncdf4::ncatt_get(nc, "time", "units")$value

  # --- read coordinates -------------------------------------------------------
  # Helper: look for a name in variables first, then in dimensions
  .read_coord <- function(nc, names) {
    for (nm in names) {
      if (nm %in% names(nc$var)) return(ncdf4::ncvar_get(nc, nm))
      if (nm %in% names(nc$dim)) return(ncdf4::ncvar_get(nc, nm))
    }
    NULL
  }

  lon  <- .read_coord(nc, c("lon", "longitude"))
  lat  <- .read_coord(nc, c("lat", "latitude"))
  rlon <- if ("rlon" %in% names(nc$dim)) ncdf4::ncvar_get(nc, "rlon") else NULL
  rlat <- if ("rlat" %in% names(nc$dim)) ncdf4::ncvar_get(nc, "rlat") else NULL

  # --- read data array --------------------------------------------------------
  dat <- ncdf4::ncvar_get(nc, var)

  # --- assemble output --------------------------------------------------------
  out <- list(
    data   = dat,
    coords = list(lat = lat, lon = lon, rlat = rlat, rlon = rlon),
    time   = list(values = time_vals, units = time_units),
    meta   = list(
      variable  = var,
      grid_type = grid_type,
      reader    = "ncdf4",
      file      = normalizePath(file, winslash = "/", mustWork = FALSE)
    )
  )

  class(out) <- c("hm_hazard", "list")
  out
}



#' Standardize spatial coordinates
#'
#' Converts the raw NetCDF coordinate vectors or matrices stored in an
#' `hm_hazard` object into a canonical data.frame with one row per grid cell.
#'
#' @details
#' Two coordinate encodings are supported:
#' \itemize{
#'   \item \strong{Regular grids:} `lat` and `lon` are 1D vectors.
#'     All combinations are expanded with [expand.grid()].
#'   \item \strong{Rotated/curvilinear grids:} `lat` and `lon` are 2D matrices
#'     (e.g., \code{lat(rlat, rlon)} as in EURO-CORDEX products).
#'     Rotated indices `rlat` and `rlon` are appended when available.
#' }
#'
#' @param x An `hm_hazard` object returned by [hm_read_netcdf()].
#'
#' @return The same `hm_hazard` object with `x$coords` replaced by a
#'   data.frame containing at least `cell_id`, `lat`, and `lon`. For rotated
#'   grids, `rlon` and `rlat` columns are also included.
#'
#' @seealso [hm_read_netcdf()], [hm_get_grid_geometry()]
#' @export
hm_standardize_coords <- function(x) {

  if (is.null(x) || !"hm_hazard" %in% class(x))
    stop("`x` must be an object of class `hm_hazard`.")
  if (is.null(x$coords) || !is.list(x$coords))
    stop("`x$coords` must be a list. Run `hm_read_netcdf()` first.")

  lat  <- x$coords$lat
  lon  <- x$coords$lon
  rlat <- x$coords$rlat
  rlon <- x$coords$rlon

  # --- regular grid: lat and lon are 1D ---------------------------------------
  is_1d <- function(z) is.null(dim(z)) || length(dim(z)) == 1L

  if (is_1d(lat) && is_1d(lon)) {
    g <- expand.grid(lat = as.vector(lat), lon = as.vector(lon))
    g$cell_id <- seq_len(nrow(g))
    x$coords  <- g[, c("cell_id", "lat", "lon")]
    return(x)
  }

  # --- rotated/curvilinear grid: lat and lon are 2D matrices -----------------
  if (is.matrix(lat) && is.matrix(lon)) {
    n <- prod(dim(lat))

    g <- data.frame(
      cell_id = seq_len(n),
      lat     = as.vector(lat),
      lon     = as.vector(lon)
    )

    if (!is.null(rlat) && !is.null(rlon)) {
      idx <- expand.grid(rlon = rlon, rlat = rlat)
      if (nrow(idx) == n) {
        g$rlon <- idx$rlon
        g$rlat <- idx$rlat
      }
    }

    x$coords <- g
    return(x)
  }

  stop("Unsupported coordinate structure in `x$coords`.")
}



#' Extract grid geometry
#'
#' Computes the grid shape, spacing, and bounding box from the coordinate
#' table stored in an `hm_hazard` object, and adds them as `x$grid`.
#'
#' @param x An `hm_hazard` object after [hm_standardize_coords()].
#'
#' @return The same `hm_hazard` object with a new `$grid` element containing:
#'   \describe{
#'     \item{`nx`, `ny`}{Number of grid cells along longitude and latitude.}
#'     \item{`dx`, `dy`}{Approximate grid spacing in degrees (`NA` for rotated).}
#'     \item{`bbox`}{Named vector: `lon_min`, `lon_max`, `lat_min`, `lat_max`.}
#'     \item{`n_cells`}{Total number of grid cells.}
#'     \item{`grid_type`}{`"regular"` or `"rotated"`.}
#'   }
#'
#' @seealso [hm_standardize_coords()]
#' @export
hm_get_grid_geometry <- function(x) {

  if (is.null(x) || !"hm_hazard" %in% class(x))
    stop("`x` must be an object of class `hm_hazard`.")
  if (is.null(x$coords) || !is.data.frame(x$coords))
    stop("`x$coords` must be a data.frame. Run `hm_standardize_coords()` first.")

  lon <- x$coords$lon
  lat <- x$coords$lat

  if (is.null(lon) || is.null(lat))
    stop("`x$coords` must contain `lon` and `lat` columns.")

  bbox <- c(
    lon_min = min(lon, na.rm = TRUE),
    lon_max = max(lon, na.rm = TRUE),
    lat_min = min(lat, na.rm = TRUE),
    lat_max = max(lat, na.rm = TRUE)
  )

  n_cells   <- nrow(x$coords)
  grid_type <- if (!is.null(x$meta$grid_type)) x$meta$grid_type else NA_character_

  # --- rotated grid: use rlat/rlon index counts, spacing is not meaningful ----
  if (!is.na(grid_type) && grid_type == "rotated") {
    if (all(c("rlon", "rlat") %in% names(x$coords))) {
      nx <- length(unique(x$coords$rlon))
      ny <- length(unique(x$coords$rlat))
    } else {
      nx <- NA_integer_
      ny <- NA_integer_
    }
    dx <- NA_real_
    dy <- NA_real_

    # --- regular grid: derive spacing from sorted unique coordinates ------------
  } else {
    lon_u <- sort(unique(lon))
    lat_u <- sort(unique(lat))

    nx <- length(lon_u)
    ny <- length(lat_u)

    dx <- if (length(lon_u) > 1) stats::median(diff(lon_u), na.rm = TRUE) else NA_real_
    dy <- if (length(lat_u) > 1) stats::median(diff(lat_u), na.rm = TRUE) else NA_real_

    # If the grid is not perfectly rectangular, dimensions are ambiguous
    if (nx * ny != n_cells) {
      nx <- NA_integer_
      ny <- NA_integer_
    }
  }

  x$grid <- list(
    nx        = nx,
    ny        = ny,
    dx        = dx,
    dy        = dy,
    bbox      = bbox,
    n_cells   = n_cells,
    grid_type = grid_type
  )

  x
}



#' Decode the NetCDF time axis
#'
#' Converts CF-style time units (e.g., `"days since 1970-01-01"`) into an R
#' `Date` or `POSIXct` vector and replaces the raw `x$time` list.
#'
#' @param x An `hm_hazard` object returned by [hm_read_netcdf()].
#' @param tz Time zone string used when units are sub-daily (default `"UTC"`).
#'
#' @return The same `hm_hazard` object with `x$time` replaced by a `Date`
#'   vector (for daily data) or a `POSIXct` vector (for sub-daily data).
#'
#' @seealso [hm_read_netcdf()], [hm_to_matrix()]
#' @export
hm_decode_time <- function(x, tz = "UTC") {

  if (is.null(x) || !"hm_hazard" %in% class(x))
    stop("`x` must be an object of class `hm_hazard`.")
  if (is.null(x$time) || !is.list(x$time) ||
      is.null(x$time$values) || is.null(x$time$units))
    stop("`x$time` must be a list with `values` and `units`. Run `hm_read_netcdf()` first.")

  vals  <- x$time$values
  units <- x$time$units

  if (!is.numeric(vals))
    stop("`x$time$values` must be numeric.")
  if (!is.character(units) || length(units) != 1L)
    stop("`x$time$units` must be a single character string.")

  u <- trimws(units)
  if (!grepl(" since ", u, fixed = TRUE))
    stop("Unsupported time units (missing ' since '): ", units)

  # --- parse "unit since origin" ----------------------------------------------
  parts     <- strsplit(u, " since ", fixed = TRUE)[[1L]]
  base_unit <- tolower(trimws(parts[[1L]]))
  origin_str <- trimws(parts[[2L]])

  # Remove trailing Z before parsing (some files use ISO 8601 format)
  origin_str <- sub("Z$", "", origin_str)
  origin     <- suppressWarnings(as.POSIXct(origin_str, tz = tz))

  if (is.na(origin))
    origin <- suppressWarnings(as.POSIXct(paste0(origin_str, " 00:00:00"), tz = tz))
  if (is.na(origin))
    stop("Could not parse time origin from units string: ", units)

  # --- convert to seconds -----------------------------------------------------
  mult <- switch(
    base_unit,
    "day"    = , "days"    = 86400,
    "hour"   = , "hours"   = 3600,
    "minute" = , "minutes" = ,
    "min"    = , "mins"    = 60,
    "second" = , "seconds" = ,
    "sec"    = , "secs"    = 1,
    NA_real_
  )

  if (!is.finite(mult))
    stop("Unsupported base time unit: `", base_unit, "`.")

  tposix <- origin + vals * mult

  # Return Date for daily data, POSIXct for sub-daily
  x$time <- if (mult == 86400) as.Date(tposix, tz = tz) else tposix
  x
}



#' Convert hazard data to a time-by-cell matrix
#'
#' Reshapes `x$data` (a multi-dimensional array) into a plain numeric matrix
#' where rows are time steps and columns are grid cells. The time vector and
#' coordinate table are attached as attributes so downstream functions do not
#' need to carry the full `hm_hazard` object.
#'
#' @param x An `hm_hazard` object after [hm_standardize_coords()] and
#'   [hm_decode_time()].
#'
#' @return A numeric matrix of dimension `T x N` (time steps by grid cells)
#'   with two attributes:
#'   \describe{
#'     \item{`time`}{The decoded time vector from `x$time`.}
#'     \item{`coords`}{The coordinate data.frame from `x$coords`.}
#'   }
#'
#' @seealso [hm_decode_time()], [hm_standardize_coords()]
#' @export
hm_to_matrix <- function(x) {

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

  dat <- x$data
  d   <- dim(dat)

  if (length(d) < 2L)
    stop("`x$data` must have at least 2 dimensions (space + time).")

  T_steps <- length(x$time)
  N_cells <- nrow(x$coords)

  # Identify which dimension corresponds to time
  tidx <- which(d == T_steps)
  if (length(tidx) != 1L)
    stop(
      "Cannot uniquely identify the time dimension. ",
      "Expected one dimension of length ", T_steps,
      ", found: ", paste(d, collapse = " x ")
    )

  # All other dimensions are spatial
  sidx <- setdiff(seq_along(d), tidx)
  if (prod(d[sidx]) != N_cells)
    stop(
      "Spatial size mismatch: product of spatial dimensions = ", prod(d[sidx]),
      " but number of coordinate rows = ", N_cells,
      ". Array dimensions: ", paste(d, collapse = " x ")
    )

  # Move time to the first dimension, then flatten to a matrix
  dat_reordered <- aperm(dat, c(tidx, sidx))
  X <- matrix(as.vector(dat_reordered), nrow = T_steps, ncol = N_cells)

  attr(X, "time")   <- x$time
  attr(X, "coords") <- x$coords
  X
}
