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




