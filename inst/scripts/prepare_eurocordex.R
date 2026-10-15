# prepare_eurocordex.R
# Merge and clip EURO-CORDEX historical wind-gust files downloaded from the
# ESGF portal into a single NetCDF ready for hm_read_netcdf().
#
# Requirements: ncdf4
#   install.packages("ncdf4")
#
# Source files follow the CORDEX naming convention:
#   wsgsmax_EUR-11_<GCM>_historical_<run>_<RCM>_<version>_day_<start>-<end>.nc
# Download from: https://esgf-data.dkrz.de

library(ncdf4)

# --- configuration -----------------------------------------------------------
HISTORICAL_DIR <- "historical"       # folder containing the raw 5-year files
OUTPUT_DIR     <- "eurocordex_prepared"

LAT_MIN <- 43.7;  LAT_MAX <- 47.2   # Northeast Italy bounding box
LON_MIN <- 6.5;   LON_MAX <- 14.2
# -----------------------------------------------------------------------------

dir.create(OUTPUT_DIR, showWarnings = FALSE)


load_and_clip <- function(data_dir) {
  files <- list.files(data_dir, pattern = "\\.nc$", full.names = TRUE)
  if (length(files) == 0L)
    stop("No NetCDF files found in: ", data_dir)

  chunks <- lapply(files, function(f) {
    nc  <- nc_open(f)
    on.exit(nc_close(nc))

    lon <- ncvar_get(nc, if ("lon" %in% names(nc$var)) "lon" else "longitude")
    lat <- ncvar_get(nc, if ("lat" %in% names(nc$var)) "lat" else "latitude")

    # CORDEX lon is sometimes [0, 360] — convert to [-180, 180]
    if (any(lon > 180)) lon <- ((lon + 180) %% 360) - 180

    in_bbox <- which(lat >= LAT_MIN & lat <= LAT_MAX &
                     lon >= LON_MIN & lon <= LON_MAX, arr.ind = TRUE)
    row_rng <- range(in_bbox[, 1])
    col_rng <- range(in_bbox[, 2])

    list(
      data   = ncvar_get(nc, "wsgsmax",
                         start = c(row_rng[1], col_rng[1], 1),
                         count = c(diff(row_rng) + 1L,
                                   diff(col_rng) + 1L, -1L)),
      time   = ncvar_get(nc, "time"),
      tunits = ncatt_get(nc, "time", "units")$value,
      lon    = lon[row_rng[1]:row_rng[2], col_rng[1]:col_rng[2]],
      lat    = lat[row_rng[1]:row_rng[2], col_rng[1]:col_rng[2]]
    )
  })

  # Concatenate along the time dimension using base R
  list(
    data   = do.call(function(...) {
               arrays <- list(...)
               dims   <- dim(arrays[[1]])
               dims[3] <- sum(vapply(arrays, function(a) dim(a)[3], integer(1)))
               out    <- array(NA_real_, dims)
               t_pos  <- 1L
               for (a in arrays) {
                 nt <- dim(a)[3]
                 out[, , t_pos:(t_pos + nt - 1L)] <- a
                 t_pos <- t_pos + nt
               }
               out
             }, lapply(chunks, `[[`, "data")),
    time   = do.call(c, lapply(chunks, `[[`, "time")),
    tunits = chunks[[1]]$tunits,
    lon    = chunks[[1]]$lon,
    lat    = chunks[[1]]$lat
  )
}


save_netcdf <- function(ds, out_file) {
  nr  <- dim(ds$lon)[1]
  nc_ <- dim(ds$lon)[2]
  nt  <- length(ds$time)

  d_rlon <- ncdim_def("rlon", "degrees", seq_len(nr))
  d_rlat <- ncdim_def("rlat", "degrees", seq_len(nc_))
  d_time <- ncdim_def("time", ds$tunits, ds$time, unlim = TRUE)

  v_lon  <- ncvar_def("lon",     "degrees_east",  list(d_rlon, d_rlat), NA)
  v_lat  <- ncvar_def("lat",     "degrees_north", list(d_rlon, d_rlat), NA)
  v_wsg  <- ncvar_def("wsgsmax", "m s-1",
                      list(d_rlon, d_rlat, d_time),
                      missval  = 1e20,
                      longname = "Daily maximum wind gust speed")

  nc_out <- nc_create(out_file, list(v_lon, v_lat, v_wsg))
  ncvar_put(nc_out, v_lon,  ds$lon)
  ncvar_put(nc_out, v_lat,  ds$lat)
  ncvar_put(nc_out, v_wsg,  ds$data)
  nc_close(nc_out)

  message("Saved: ", basename(out_file))
}


# Merge and clip the historical files (1970-2005)
message("Processing historical files in: ", HISTORICAL_DIR)
ds <- load_and_clip(HISTORICAL_DIR)

out_file <- file.path(OUTPUT_DIR, "EuroCordex_wsgsmax_1970_2005.nc")
save_netcdf(ds, out_file)

message("Ready for: hm_read_netcdf('", out_file, "', var = 'wsgsmax')")
