#' Read a hazard NetCDF dataset
#'
#' Reads a NetCDF file containing a hazard variable on a spatial grid
#' (regular or rotated) and returns an `hm_hazard` object. CRS is not
#' required.
#'
#' @details
#' Supports CF-convention NetCDF files with regular grids (e.g. Copernicus,
#' ERA5) and rotated or curvilinear grids (e.g. EURO-CORDEX). Coordinate
#' names \code{lat}/\code{lon} and \code{latitude}/\code{longitude} are
#' detected automatically.
#'
#' @param file Path to a NetCDF file.
#' @param var Name of the hazard variable. If \code{NULL} and the file
#'   contains exactly one variable, it is selected automatically.
#'
#' @return An object of class \code{hm_hazard}.
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

  var <- .hm_resolve_var(nc, var)

  dnames    <- names(nc$dim)
  grid_type <- if (all(c("rlat", "rlon") %in% dnames)) "rotated" else "regular"

  if (!"time" %in% dnames)
    stop("No `time` dimension found in the NetCDF file.")

  time_vals  <- ncdf4::ncvar_get(nc, "time")
  time_units <- ncdf4::ncatt_get(nc, "time", "units")$value

  lon  <- .hm_read_coord(nc, c("lon", "longitude"))
  lat  <- .hm_read_coord(nc, c("lat", "latitude"))
  rlon <- if ("rlon" %in% dnames) ncdf4::ncvar_get(nc, "rlon") else NULL
  rlat <- if ("rlat" %in% dnames) ncdf4::ncvar_get(nc, "rlat") else NULL

  out <- list(
    data   = ncdf4::ncvar_get(nc, var),
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

.hm_resolve_var <- function(nc, var) {
  vnames <- names(nc$var)
  if (is.null(var)) {
    if (length(vnames) != 1L)
      stop("Multiple variables found. Specify `var`. Candidates: ",
           paste(vnames, collapse = ", "))
    return(vnames[[1L]])
  }
  if (!var %in% vnames)
    stop("Variable `", var, "` not found. Candidates: ",
         paste(vnames, collapse = ", "))
  var
}

.hm_read_coord <- function(nc, names) {
  for (nm in names) {
    if (nm %in% names(nc$var)) return(ncdf4::ncvar_get(nc, nm))
    if (nm %in% names(nc$dim)) return(ncdf4::ncvar_get(nc, nm))
  }
  NULL
}



#' Standardize spatial coordinates
#'
#' Converts raw NetCDF coordinates into a canonical data.frame with one row
#' per grid cell.
#'
#' @details
#' Regular grids supply \code{lat} and \code{lon} as 1D vectors; all
#' combinations are expanded with \code{expand.grid()}. Rotated or
#' curvilinear grids supply \code{lat} and \code{lon} as 2D matrices, with
#' rotated indices appended when available.
#'
#' @param x An \code{hm_hazard} object returned by [hm_read_netcdf()].
#'
#' @return The same object with \code{x$coords} replaced by a data.frame
#'   containing \code{cell_id}, \code{lat}, and \code{lon}.
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

  is_1d <- function(z) is.null(dim(z)) || length(dim(z)) == 1L

  if (is_1d(lat) && is_1d(lon)) {
    g         <- expand.grid(lat = as.vector(lat), lon = as.vector(lon))
    g$cell_id <- seq_len(nrow(g))
    x$coords  <- g[, c("cell_id", "lat", "lon")]
    return(x)
  }

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
#' Computes the bounding box, grid spacing, and cell count from the
#' coordinate table and stores them in \code{x$grid}.
#'
#' @param x An \code{hm_hazard} object after [hm_standardize_coords()].
#'
#' @return The same object with a new \code{$grid} element containing
#'   \code{nx}, \code{ny}, \code{dx}, \code{dy}, \code{bbox},
#'   \code{n_cells}, and \code{grid_type}.
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

  bbox <- c(
    lon_min = min(lon, na.rm = TRUE),
    lon_max = max(lon, na.rm = TRUE),
    lat_min = min(lat, na.rm = TRUE),
    lat_max = max(lat, na.rm = TRUE)
  )

  n_cells   <- nrow(x$coords)
  grid_type <- x$meta$grid_type %||% NA_character_

  if (!is.na(grid_type) && grid_type == "rotated") {
    nx <- length(unique(x$coords$rlon))
    ny <- length(unique(x$coords$rlat))
    dx <- NA_real_
    dy <- NA_real_
  } else {
    lon_u <- sort(unique(lon))
    lat_u <- sort(unique(lat))
    nx    <- length(lon_u)
    ny    <- length(lat_u)
    dx    <- if (nx > 1) stats::median(diff(lon_u)) else NA_real_
    dy    <- if (ny > 1) stats::median(diff(lat_u)) else NA_real_
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

`%||%` <- function(x, y) if (!is.null(x)) x else y



#' Decode the NetCDF time axis
#'
#' Converts CF-style time values (e.g. \code{"days since 1970-01-01"}) to
#' an R \code{Date} or \code{POSIXct} vector.
#'
#' @param x An \code{hm_hazard} object returned by [hm_read_netcdf()].
#' @param tz Time zone string for sub-daily data. Defaults to \code{"UTC"}.
#'
#' @return The same object with \code{x$time} replaced by a \code{Date}
#'   vector (daily data) or \code{POSIXct} vector (sub-daily data).
#'
#' @seealso [hm_read_netcdf()], [hm_to_matrix()]
#' @export
hm_decode_time <- function(x, tz = "UTC") {

  if (is.null(x) || !"hm_hazard" %in% class(x))
    stop("`x` must be an object of class `hm_hazard`.")
  if (!is.list(x$time) || is.null(x$time$values) || is.null(x$time$units))
    stop("`x$time` must be a list with `values` and `units`. ",
         "Run `hm_read_netcdf()` first.")

  vals  <- x$time$values
  units <- x$time$units

  if (!grepl(" since ", units, fixed = TRUE))
    stop("Unsupported time units (missing ' since '): ", units)

  parts     <- strsplit(trimws(units), " since ", fixed = TRUE)[[1L]]
  base_unit <- tolower(trimws(parts[[1L]]))
  origin    <- .hm_parse_origin(trimws(parts[[2L]]), tz)

  mult <- .hm_time_multiplier(base_unit)
  if (!is.finite(mult))
    stop("Unsupported base time unit: `", base_unit, "`.")

  tposix   <- origin + vals * mult
  x$time   <- if (mult == 86400) as.Date(tposix, tz = tz) else tposix
  x
}

.hm_parse_origin <- function(origin_str, tz) {
  origin_str <- sub("Z$", "", origin_str)
  origin     <- suppressWarnings(as.POSIXct(origin_str, tz = tz))
  if (is.na(origin))
    origin <- suppressWarnings(
      as.POSIXct(paste0(origin_str, " 00:00:00"), tz = tz)
    )
  if (is.na(origin))
    stop("Could not parse time origin: ", origin_str)
  origin
}

.hm_time_multiplier <- function(base_unit) {
  switch(base_unit,
    "day"    = , "days"    = 86400,
    "hour"   = , "hours"   = 3600,
    "minute" = , "minutes" = ,
    "min"    = , "mins"    = 60,
    "second" = , "seconds" = ,
    "sec"    = , "secs"    = 1,
    NA_real_
  )
}



#' Convert hazard array to a time-by-cell matrix
#'
#' Reshapes \code{x$data} into a \eqn{T \times N} numeric matrix where rows
#' are time steps and columns are grid cells.
#'
#' @param x An \code{hm_hazard} object after [hm_standardize_coords()] and
#'   [hm_decode_time()].
#'
#' @return A numeric matrix with attributes \code{time} and \code{coords}.
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
  if (!inherits(x$time, c("Date", "POSIXct", "POSIXt")))
    stop("`x$time` must be decoded. Run `hm_decode_time()` first.")

  d       <- dim(x$data)
  T_steps <- length(x$time)
  N_cells <- nrow(x$coords)

  tidx <- which(d == T_steps)
  if (length(tidx) != 1L)
    stop("Cannot identify the time dimension uniquely. ",
         "Expected one dimension of length ", T_steps,
         ", found: ", paste(d, collapse = " x "))

  sidx <- setdiff(seq_along(d), tidx)
  if (prod(d[sidx]) != N_cells)
    stop("Spatial size mismatch: ", prod(d[sidx]),
         " cells in data vs ", N_cells, " rows in coords.")

  X <- matrix(
    as.vector(aperm(x$data, c(tidx, sidx))),
    nrow = T_steps,
    ncol = N_cells
  )

  attr(X, "time")   <- x$time
  attr(X, "coords") <- x$coords
  X
}



#' Remove sea and no-data cells from a hazard matrix
#'
#' Drops columns where the fraction of finite positive values falls below
#' \code{min_valid}. Needed for rotated-grid datasets such as EURO-CORDEX
#' where sea cells produce all-NA columns.
#'
#' @param X A numeric matrix as returned by [hm_to_matrix()].
#' @param min_valid Minimum fraction of valid observations required to keep
#'   a column. Defaults to \code{0.8}.
#'
#' @return A matrix with only land columns retained. Attributes
#'   \code{land_idx}, \code{coords}, and \code{time} are attached.
#'
#' @seealso [hm_to_matrix()], [hm_fit_marginals()]
#' @export
hm_filter_land_cells <- function(X, min_valid = 0.8) {

  if (!is.matrix(X) || !is.numeric(X))
    stop("`X` must be a numeric matrix.")
  if (!is.numeric(min_valid) || length(min_valid) != 1L ||
      min_valid <= 0 || min_valid > 1)
    stop("`min_valid` must be a single numeric value in (0, 1].")

  valid_frac <- colSums(is.finite(X) & X > 0) / nrow(X)
  land_idx   <- which(valid_frac >= min_valid)

  if (length(land_idx) == 0L)
    stop("No columns meet `min_valid` = ", min_valid,
         ". Consider lowering the threshold.")

  n_dropped <- ncol(X) - length(land_idx)
  if (n_dropped > 0L)
    message(n_dropped, " column(s) removed. ",
            length(land_idx), " land cells retained.")

  X_land <- X[, land_idx, drop = FALSE]

  attr(X_land, "land_idx") <- land_idx
  attr(X_land, "time")     <- attr(X, "time")

  coords <- attr(X, "coords")
  if (!is.null(coords) && is.data.frame(coords))
    attr(X_land, "coords") <- coords[land_idx, , drop = FALSE]

  X_land
}
