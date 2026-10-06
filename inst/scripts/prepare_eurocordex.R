# prepare_eurocordex.R
# Merge and clip EURO-CORDEX 5-year NetCDF files downloaded from the ESGF portal.
# Output: NetCDF files ready for hm_read_netcdf(file, var = "wsgsmax")
#
# Requirements: ncdf4
#   install.packages("ncdf4")
#
# Source files follow the CORDEX naming convention:
#   wsgsmax_EUR-11_<GCM>_<scenario>_<run>_<RCM>_<version>_day_<start>-<end>.nc
# Download from: https://esgf-data.dkrz.de

library(ncdf4)

# --- configuration -----------------------------------------------------------
HISTORICAL_DIR <- "historical"
RCP45_DIR      <- "rcp4_5"
RCP85_DIR      <- "rcp8_5"
OUTPUT_DIR     <- "eurocordex_prepared"

LAT_MIN <- 43.7;  LAT_MAX <- 47.2   # Northeast Italy
LON_MIN <- 6.5;   LON_MAX <- 14.2

PERIODS <- list(
  near_future = c("2021-01-01", "2040-12-31"),
  mid_century = c("2041-01-01", "2070-12-31"),
  far_future  = c("2071-01-01", "2100-12-31")
)
# -----------------------------------------------------------------------------

dir.create(OUTPUT_DIR, showWarnings = FALSE)


load_and_clip <- function(data_dir) {
  files <- list.files(data_dir, pattern = "\\.nc$", full.names = TRUE)
  if (length(files) == 0L)
    stop("No NetCDF files found in: ", data_dir)

  # Read and concatenate along time
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

  list(
    data   = abind::abind(lapply(chunks, `[[`, "data"), along = 3),
    time   = do.call(c, lapply(chunks, `[[`, "time")),
    tunits = chunks[[1]]$tunits,
    lon    = chunks[[1]]$lon,
    lat    = chunks[[1]]$lat
  )
}


save_period <- function(ds, start_date, end_date, out_file) {
  # Convert time to Date for subsetting
  origin <- as.Date(sub(".*since ", "", ds$tunits))
  unit   <- trimws(sub(" since.*", "", ds$tunits))
  mult   <- switch(unit, days = 1, hours = 1/24, seconds = 1/86400, 1)
  dates  <- origin + ds$time * mult

  keep <- which(dates >= as.Date(start_date) & dates <= as.Date(end_date))
  if (length(keep) == 0L) {
    message("No data in period ", start_date, " — ", end_date, ". Skipping.")
    return(invisible(NULL))
  }

  nr <- dim(ds$lon)[1]
  nc_ <- dim(ds$lon)[2]

  d_rlon <- ncdim_def("rlon", "degrees", seq_len(nr))
  d_rlat <- ncdim_def("rlat", "degrees", seq_len(nc_))
  d_time <- ncdim_def("time", ds$tunits, ds$time[keep], unlim = TRUE)

  v_lon  <- ncvar_def("lon", "degrees_east",  list(d_rlon, d_rlat), NA)
  v_lat  <- ncvar_def("lat", "degrees_north", list(d_rlon, d_rlat), NA)
  v_wsg  <- ncvar_def("wsgsmax", "m s-1",
                      list(d_rlon, d_rlat, d_time),
                      missval = 1e20,
                      longname = "Daily maximum wind gust speed")

  nc_out <- nc_create(out_file, list(v_lon, v_lat, v_wsg))
  ncvar_put(nc_out, v_lon,  ds$lon)
  ncvar_put(nc_out, v_lat,  ds$lat)
  ncvar_put(nc_out, v_wsg,  ds$data[, , keep])
  nc_close(nc_out)

  message("Saved: ", basename(out_file))
}


# Historical baseline (1970–2005)
message("Processing historical ...")
ds_hist <- load_and_clip(HISTORICAL_DIR)
save_period(ds_hist, "1970-01-01", "2005-12-31",
            file.path(OUTPUT_DIR, "EuroCordex_wsgsmax_1970_2005.nc"))
message("Ready for: hm_read_netcdf('EuroCordex_wsgsmax_1970_2005.nc', var = 'wsgsmax')")

# Climate projections
for (scenario in c("RCP45", "RCP85")) {
  data_dir <- if (scenario == "RCP45") RCP45_DIR else RCP85_DIR
  message("\nProcessing ", scenario, " ...")
  ds_proj <- load_and_clip(data_dir)

  for (period in names(PERIODS)) {
    save_period(ds_proj, PERIODS[[period]][1], PERIODS[[period]][2],
                file.path(OUTPUT_DIR, paste0(scenario, "_", period, ".nc")))
  }
}
