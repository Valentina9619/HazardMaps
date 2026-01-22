#' Read a hazard NetCDF dataset
#'
#' Reads a NetCDF file containing a hazard variable on a spatial grid (regular or rotated)
#' and returns a basic `hm_hazard` object. CRS is NOT required.
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

                                              lon <- if ("lon" %in% names(nc$var)) ncdf4::ncvar_get(nc, "lon") else if ("longitude" %in% names(nc$var)) ncdf4::ncvar_get(nc, "longitude") else NULL
                                              lat <- if ("lat" %in% names(nc$var)) ncdf4::ncvar_get(nc, "lat") else if ("latitude" %in% names(nc$var)) ncdf4::ncvar_get(nc, "latitude") else NULL

                                              rlon <- if ("rlon" %in% names(nc$dim)) ncdf4::ncvar_get(nc, "rlon") else NULL
                                              rlat <- if ("rlat" %in% names(nc$dim)) ncdf4::ncvar_get(nc, "rlat") else NULL

                                              dat <- ncdf4::ncvar_get(nc, var)

                                              out <- list(data = dat,
                                                          coords = list(lat = lat, lon = lon, rlat = rlat, rlon = rlon),
                                                          time = list(values = time_vals, units = time_units),
                                                          meta = list(variable = var,
                                                          grid_type = grid_type,
                                                          reader = "ncdf4",
                                                          file = normalizePath(file, winslash = "/", mustWork = FALSE)))
                                              class(out) <- c("hm_hazard", "list")
                                              out}



#' Standardize spatial coordinates
#'
#' Converts raw NetCDF coordinates (regular 1D lat/lon or rotated/curvilinear 2D lat/lon)
#' into a canonical coordinate table with one row per grid cell.
#'
#' @param x An `hm_hazard` object returned by [hm_read_netcdf()].
#'
#' @return The same `hm_hazard` object, but with `x$coords` replaced by a data.frame
#'   containing at least `cell_id`, `lat`, and `lon`. If available, `rlat` and `rlon`
#'   are also included.
#' @export
hm_standardize_coords <- function(x) {if (is.null(x) || !"hm_hazard" %in% class(x)) stop("`x` must be an object of class `hm_hazard`.")
                                      if (is.null(x$coords) || !is.list(x$coords)) stop("`x$coords` must be a list with lat/lon.")
                                      lat <- x$coords$lat; lon <- x$coords$lon
                                      if (is.null(lat) || is.null(lon)) stop("Missing lat/lon in `x$coords`. Supported names are lat/lon or latitude/longitude in the NetCDF.")
                                      if (!is.numeric(lat) || !is.numeric(lon)) stop("lat/lon must be numeric.")

                                      rlat <- x$coords$rlat; rlon <- x$coords$rlon

                                      # Case 1: regular grid (lat/lon 1D)
                                      if (is.vector(lat) && is.vector(lon) && length(dim(lat)) == 0L && length(dim(lon)) == 0L) {g <- expand.grid(lat = lat, lon = lon)
                                                                                                                                 g$cell_id <- seq_len(nrow(g))
                                                                                                                                 g <- g[, c("cell_id","lat","lon")]
                                                                                                                                 x$coords <- g
                                                                                                                                 return(x)}

                                      # Case 2: rotated/curvilinear grid (lat/lon 2D arrays)
                                      if (length(dim(lat)) == 2L && length(dim(lon)) == 2L) {dlat <- dim(lat); dlon <- dim(lon)
                                                                                             if (any(dlat != dlon)) stop("lat and lon dimensions do not match.")
                                                                                             n <- dlat[[1L]] * dlat[[2L]]
                                                                                             g <- data.frame(cell_id = seq_len(n),
                                                                                                             lat = as.vector(lat),
                                                                                                             lon = as.vector(lon))
                                                                                             # attach rlat/rlon if present and consistent
                                                                                             if (!is.null(rlat) && !is.null(rlon) && is.vector(rlat) && is.vector(rlon) && length(rlat) == dlat[[1L]] && length(rlon) == dlat[[2L]]) {idx <- expand.grid(rlat = rlat, rlon = rlon)
                                                                                                                                                                                                                                      g$rlat <- idx$rlat
                                                                                                                                                                                                                                      g$rlon <- idx$rlon}
                                                                                              x$coords <- g
                                                                                              return(x)}

                                      stop("Unsupported lat/lon structure. Expected 1D (regular grid) or 2D (rotated/curvilinear grid).")}



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



#' Convert hazard array to a time X cell matrix
#'
#' Converts `x$data` (an array with one time dimension) into a numeric matrix with:
#' rows = time steps, columns = grid cells (cell_id). Time and coordinates are stored as attributes.
#'
#' @param x An `hm_hazard` object after [hm_standardize_coords()] and [hm_decode_time()].
#'
#' @return A numeric matrix of dimension (T X N) with attributes:
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







