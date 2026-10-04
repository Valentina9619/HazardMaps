#' Read a hazard NetCDF dataset
#'
#' Reads a NetCDF file containing a hazard variable on a spatial grid (regular or rotated)
#' and returns a basic `hm_hazard` object. CRS is NOT required.
#'
#' @details
#' The function supports NetCDF files following common CF conventions for spatial
#' coordinates. Latitude and longitude may be provided either as coordinate variables
#' or as dimension coordinates, and may be encoded as:
#' \itemize{
#'   \item one-dimensional regular grids (e.g., Copernicus / ERA5 products),
#'   \item two-dimensional curvilinear or rotated grids (e.g., EURO-CORDEX products).
#' }
#' Coordinate names \code{lat}/\code{lon} and \code{latitude}/\code{longitude}
#' are automatically detected.
#'
#' @param file Path to a NetCDF file.
#' @param var Optional. Name of the hazard variable. If `NULL` and multiple variables exist,
#'   an error is raised (safe default).
#'
#' @return An object of class `hm_hazard`.
#' @export
hm_read_netcdf <- function(file, var = NULL) {if (!is.character(file) || length(file) != 1L) stop("`file` must be a single path string.")
                                              if (!file.exists(file)) stop("File not found: ", file)

                                              nc <- ncdf4::nc_open(file); on.exit(ncdf4::nc_close(nc), add = TRUE)

                                              vnames <- names(nc$var)
                                              if (is.null(var)) {if (length(vnames) != 1L) stop("Multiple variables found. Specify `var`. Candidates: ", paste(vnames, collapse = ", "))
                                                                 var <- vnames[[1L]]
                                                                  } else {if (!var %in% vnames) stop("Variable `", var, "` not found. Candidates: ", paste(vnames, collapse = ", "))}

                                              dnames <- names(nc$dim)
                                              if (!("time" %in% dnames)) stop("No `time` dimension found in NetCDF.")
                                              grid_type <- if (all(c("rlat","rlon") %in% dnames)) "rotated" else "regular"

                                              time_vals <- ncdf4::ncvar_get(nc, "time")
                                              time_units <- ncdf4::ncatt_get(nc, "time", "units")$value


                                              lon <- if ("lon" %in% names(nc$var)) {ncdf4::ncvar_get(nc, "lon")} else if ("longitude" %in% names(nc$var)) {ncdf4::ncvar_get(nc, "longitude")}
                                                else if ("lon" %in% names(nc$dim)) {ncdf4::ncvar_get(nc, "lon")} else if ("longitude" %in% names(nc$dim)) {ncdf4::ncvar_get(nc, "longitude")}
                                                else {NULL}

                                              lat <- if ("lat" %in% names(nc$var)) {ncdf4::ncvar_get(nc, "lat")} else if ("latitude" %in% names(nc$var)) {ncdf4::ncvar_get(nc, "latitude")}
                                                else if ("lat" %in% names(nc$dim)) {ncdf4::ncvar_get(nc, "lat")} else if ("latitude" %in% names(nc$dim)) {ncdf4::ncvar_get(nc, "latitude")}
                                                else {NULL}

                                              rlon <- if ("rlon" %in% names(nc$dim)) ncdf4::ncvar_get(nc, "rlon") else NULL
                                              rlat <- if ("rlat" %in% names(nc$dim)) ncdf4::ncvar_get(nc, "rlat") else NULL

                                              dat <- ncdf4::ncvar_get(nc, var)

                                              out <- list(data = dat,
                                                          coords = list(lat = lat, lon = lon, rlat = rlat, rlon = rlon),
                                                          time = list(values = time_vals, units = time_units),
                                                          meta = list(variable = var,
                                                                      grid_type = grid_type,
                                                                      reader = "ncdf4",
                                                                      file = normalizePath(file, winslash = "/", mustWork = FALSE))
                                                          )
                                              class(out) <- c("hm_hazard", "list")
                                              out}



#' Standardize spatial coordinates
#'
#' Converts raw NetCDF coordinates into a canonical coordinate table with one row
#' per grid cell.
#'
#' @details
#' Two coordinate encodings are supported:
#' \itemize{
#'   \item \strong{Regular grids:} latitude and longitude provided as 1D coordinates
#'   (either pure vectors or 1D arrays with dimension attributes).
#'   \item \strong{Curvilinear/rotated grids:} latitude and longitude provided as 2D arrays
#'   (e.g., \code{lat(rlat, rlon)} and \code{lon(rlat, rlon)}).
#' }
#'
#' @param x An `hm_hazard` object returned by [hm_read_netcdf()].
#'
#' @return The same `hm_hazard` object, but with `x$coords` replaced by a data.frame
#'   containing at least `cell_id`, `lat`, and `lon`. Additional columns may be included
#'   when available (e.g., rotated-grid indices).
#' @export
hm_standardize_coords <- function(x) {if (is.null(x) || !"hm_hazard" %in% class(x)) stop("`x` must be an object of class `hm_hazard`.")
                                      if (is.null(x$coords) || !is.list(x$coords)) stop("`x$coords` must be a list with lat/lon.")

                                      lat <- x$coords$lat; lon <- x$coords$lon
                                      rlat <- x$coords$rlat; rlon <- x$coords$rlon

                                      if ((is.null(dim(lat)) || length(dim(lat)) == 1) && (is.null(dim(lon)) || length(dim(lon)) == 1)) {latv <- as.vector(lat)
                                                                                                                                         lonv <- as.vector(lon)

                                                                                                                                         g <- expand.grid(lat = latv, lon = lonv)
                                                                                                                                         g$cell_id <- seq_len(nrow(g))
                                                                                                                                         g <- g[, c("cell_id","lat","lon")]

                                                                                                                                         x$coords <- g
                                                                                                                                         return(x)}

                                       if (is.matrix(lat) && is.matrix(lon)) {d <- dim(lat)
                                                                              n <- d[1] * d[2]

                                                                              g <- data.frame(cell_id = seq_len(n),
                                                                                              lat = as.vector(lat),
                                                                                              lon = as.vector(lon))

                                                                              if (!is.null(rlat) && !is.null(rlon)) {idx <- expand.grid(rlon = rlon, rlat = rlat)
                                                                                                                     if (nrow(idx) == n) {g$rlon <- idx$rlon
                                                                                                                                          g$rlat <- idx$rlat}
                                                                                                                    }

                                                                               x$coords <- g
                                                                               return(x)}

                                       stop("Unsupported coordinate structure.")}



#' Extract grid geometry from hazard object
#'
#' Computes spatial grid information required for visualization.
#' This function derives the grid shape, spacing and bounding box
#' from the coordinate information contained in the hazard object.
#'
#' @param x An `hm_hazard` object after [hm_standardize_coords()].
#'
#' @return The input object with an additional `grid` element
#' containing:
#'   - `nx`, `ny` grid dimensions
#'   - `dx`, `dy` approximate grid spacing (NA for rotated grids)
#'   - `bbox` bounding box of the spatial domain
#'   - `n_cells` total number of grid cells
#'   - `grid_type` grid classification ("regular" or "rotated")
#'
#' @export
hm_get_grid_geometry <- function(x) {if (is.null(x) || !"hm_hazard" %in% class(x))
                                      stop("`x` must be an object of class `hm_hazard`.")

                                     if (is.null(x$coords) || !is.data.frame(x$coords))
                                      stop("`x$coords` must be a data.frame. Run `hm_standardize_coords()` first.")

                                     lon <- x$coords$lon
                                     lat <- x$coords$lat

                                     if (is.null(lon) || is.null(lat))
                                      stop("Coordinates must contain `lon` and `lat`.")

                                     bbox <- c(lon_min = min(lon, na.rm = TRUE),
                                               lon_max = max(lon, na.rm = TRUE),
                                               lat_min = min(lat, na.rm = TRUE),
                                               lat_max = max(lat, na.rm = TRUE))

                                     n_cells <- nrow(x$coords)

                                     grid_type <- if (!is.null(x$meta$grid_type)) x$meta$grid_type else NA_character_

                                     if (!is.na(grid_type) && grid_type == "rotated")
                                       {if (all(c("rlon","rlat") %in% names(x$coords))) {nx <- length(unique(x$coords$rlon))
                                                                                        ny <- length(unique(x$coords$rlat))
                                                                                        } else {nx <- NA_integer_
                                                                                                ny <- NA_integer_}

                                        dlon <- NA_real_
                                        dlat <- NA_real_
                                     }
                                     else {lon_u <- sort(unique(lon))
                                           lat_u <- sort(unique(lat))

                                           nx <- length(lon_u)
                                           ny <- length(lat_u)

                                           dlon <- if (length(lon_u) > 1) stats::median(diff(lon_u), na.rm = TRUE) else NA_real_
                                           dlat <- if (length(lat_u) > 1) stats::median(diff(lat_u), na.rm = TRUE) else NA_real_

                                            if (nx * ny != n_cells) {nx <- NA_integer_
                                                                     ny <- NA_integer_}
                                          }

                                     x$grid <- list(nx = nx, ny = ny,
                                                    dx = dlon, dy = dlat,
                                                    bbox = bbox,
                                                    n_cells = n_cells,
                                                    grid_type = grid_type)

                                     x}



#' Decode NetCDF time axis
#'
#' Converts CF-style NetCDF time units (e.g., "days since 1970-01-01") into
#' an R Date (or POSIXct) vector and stores it in `x$time`.
#'
#' @param x An `hm_hazard` object returned by [hm_read_netcdf()].
#' @param tz Timezone used if units are in hours/minutes/seconds (default "UTC").
#'
#' @return The same `hm_hazard` object with `x$time` replaced by a Date or POSIXct vector.
#' @export
hm_decode_time <- function(x, tz = "UTC") {if (is.null(x) || !"hm_hazard" %in% class(x)) stop("`x` must be an object of class `hm_hazard`.")
                                           if (is.null(x$time) || !is.list(x$time) || is.null(x$time$values) || is.null(x$time$units)) stop("`x$time` must be a list with `values` and `units`.")
                                           vals <- x$time$values; units <- x$time$units
                                           if (!is.numeric(vals)) stop("`x$time$values` must be numeric.")
                                           if (!is.character(units) || length(units) != 1L) stop("`x$time$units` must be a single string.")
                                           u <- trimws(units)
                                           if (!grepl(" since ", u, fixed = TRUE)) stop("Unsupported time units format (missing ' since '): ", units)

                                           parts <- strsplit(u, " since ", fixed = TRUE)[[1L]]
                                           base_unit <- tolower(trimws(parts[[1L]]))
                                           origin_str <- trimws(parts[[2L]])

                                           origin_str <- sub("Z$", "", origin_str)
                                           origin <- suppressWarnings(as.POSIXct(origin_str, tz = tz))
                                           if (is.na(origin)) origin <- suppressWarnings(as.POSIXct(paste0(origin_str, " 00:00:00"), tz = tz))
                                           if (is.na(origin)) stop("Could not parse time origin from units: ", units)

                                           mult <- if (base_unit %in% c("days", "day")) 86400 else if (base_unit %in% c("hours", "hour")) 3600 else if (base_unit %in% c("minutes", "minute", "mins", "min")) 60 else if (base_unit %in% c("seconds", "second", "secs", "sec")) 1 else NA_real_
                                           if (!is.finite(mult)) stop("Unsupported base time unit: ", base_unit)

                                           tposix <- origin + vals * mult
                                           if (mult == 86400) x$time <- as.Date(tposix, tz = tz) else x$time <- tposix
                                           x}



#' Convert hazard array to a time x cell matrix
#'
#' Converts `x$data` (an array with one time dimension) into a numeric matrix with:
#' rows = time steps, columns = grid cells (cell_id). Time and coordinates are stored as attributes.
#'
#' @param x An `hm_hazard` object after [hm_standardize_coords()] and [hm_decode_time()].
#'
#' @return A numeric matrix of dimension (T x N) with attributes:
#'   - `time`: the decoded time vector
#'   - `coords`: the coordinate table (must include `cell_id`, `lat`, `lon`)
#' @export
hm_to_matrix <- function(x) {if (is.null(x) || !"hm_hazard" %in% class(x)) stop("`x` must be an object of class `hm_hazard`.")
                             if (is.null(x$data)) stop("`x$data` is missing.")
                             if (is.null(x$coords) || !is.data.frame(x$coords)) stop("`x$coords` must be a data.frame. Run `hm_standardize_coords()` first.")
                             if (is.null(x$time) || !(inherits(x$time, "Date") || inherits(x$time, "POSIXct") || inherits(x$time, "POSIXt"))) stop("`x$time` must be decoded (Date/POSIXct). Run `hm_decode_time()` first.")

                             dat <- x$data
                             if (!is.array(dat)) stop("`x$data` must be an array.")
                             d <- dim(dat); if (is.null(d) || length(d) < 2L) stop("`x$data` must have at least 2 dimensions (time + space).")

                             T <- length(x$time); N <- nrow(x$coords)
                             tidx <- which(d == T)
                             if (length(tidx) != 1L) stop("Could not uniquely identify the time dimension in `x$data`. Expected exactly one dimension of length ", T, ". Found: ", paste(d, collapse = " x "))

                             sidx <- setdiff(seq_along(d), tidx)
                             if (prod(d[sidx]) != N) stop("Spatial size mismatch: prod(spatial dims) = ", prod(d[sidx]), " but nrow(coords) = ", N, ". dims(data) = ", paste(d, collapse = " x "))

                             dat2 <- aperm(dat, c(tidx, sidx))
                             X <- matrix(as.vector(dat2), nrow = T, ncol = N, byrow = FALSE)

                             attr(X, "time") <- x$time
                             attr(X, "coords") <- x$coords
                             X}



#' Remove sea and no-data cells from a hazard matrix
#'
#' Drops columns from the hazard matrix that contain too few valid
#' (finite, positive) observations to support marginal fitting or spatial
#' dependence estimation. This is the standard pre-processing step before
#' [hm_fit_marginals()] or [hm_empirical_kendall()] on rotated-grid datasets
#' such as EURO-CORDEX, where sea cells produce all-NA or all-zero columns.
#'
#' @details
#' A column is retained when the fraction of finite, strictly positive
#' values is at least `min_valid`. The default of `0.8` means a cell must
#' have valid data on at least 80 percent of the extreme-event days to be
#' kept. The indices of retained columns are stored as attribute
#' `"land_idx"` on the output matrix, so results can be mapped back to the
#' full coordinate table for plotting.
#'
#' @param X A numeric matrix of dimension \eqn{T \times N} as returned by
#'   [hm_to_matrix()].
#' @param min_valid A single numeric value in \eqn{(0, 1]}. Columns where
#'   the fraction of finite positive values is below this threshold are
#'   dropped. Defaults to `0.8`.
#'
#' @return A numeric matrix with only land columns retained. Two attributes
#'   are attached:
#'   \describe{
#'     \item{`land_idx`}{Integer vector of the retained column indices in
#'       the original matrix.}
#'     \item{`coords`}{Coordinate data.frame restricted to the retained
#'       columns, when the original matrix carries a `coords` attribute.}
#'   }
#'
#' @seealso [hm_to_matrix()], [hm_fit_marginals()],
#'   [hm_empirical_kendall()]
#'
#' @examples
#' \dontrun{
#' f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
#'                  package = "HazardMaps")
#' x <- hm_read_netcdf(f, var = "wsgsmax")
#' x <- hm_standardize_coords(x)
#' x <- hm_decode_time(x)
#' x <- hm_select_extreme_events(x, threshold = 25)
#' X <- hm_to_matrix(x)
#'
#' X_land <- hm_filter_land_cells(X)
#' cat("Cells before:", ncol(X), " — after:", ncol(X_land), "\n")
#' attr(X_land, "land_idx")
#' }
#'
#' @export
hm_filter_land_cells <- function(X, min_valid = 0.8) {

  if (!is.matrix(X) || !is.numeric(X))
    stop("`X` must be a numeric matrix.")
  if (!is.numeric(min_valid) || length(min_valid) != 1L ||
      min_valid <= 0 || min_valid > 1)
    stop("`min_valid` must be a single numeric value in (0, 1].")

  T_steps <- nrow(X)

  # Fraction of finite positive values per column
  valid_frac <- colSums(is.finite(X) & X > 0) / T_steps

  land_idx <- which(valid_frac >= min_valid)

  if (length(land_idx) == 0L)
    stop(
      "No columns meet the `min_valid` threshold of ", min_valid, ". ",
      "Consider lowering the threshold."
    )

  n_dropped <- ncol(X) - length(land_idx)
  if (n_dropped > 0L)
    message(n_dropped, " sea/no-data column(s) removed. ",
            length(land_idx), " land cells retained.")

  X_land <- X[, land_idx, drop = FALSE]

  attr(X_land, "land_idx") <- land_idx

  # Carry forward the coords attribute restricted to land cells
  coords <- attr(X, "coords")
  if (!is.null(coords) && is.data.frame(coords))
    attr(X_land, "coords") <- coords[land_idx, , drop = FALSE]

  # Carry forward the time attribute unchanged
  attr(X_land, "time") <- attr(X, "time")

  X_land
}
