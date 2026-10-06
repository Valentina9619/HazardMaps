# prepare_copernicus.R
# Download and merge ERA5 daily-maximum wind-gust files from Copernicus CDS.
# Output: one NetCDF ready for hm_read_netcdf(file, var = "i10fg")
#
# Requirements: ecmwfr, ncdf4
#   install.packages(c("ecmwfr", "ncdf4"))
#
# Before running: register at https://cds.climate.copernicus.eu and store
#   your credentials once with ecmwfr::wf_set_key()

library(ecmwfr)
library(ncdf4)

# --- configuration -----------------------------------------------------------
YEARS       <- 1994:2021
AREA        <- c(47.25, 10.75, 43.75, 14.25)  # N / W / S / E
DOWNLOAD_DIR <- "era5_monthly"
OUTPUT_FILE  <- "Copernicus_dailymax_1994_2021.nc"
# -----------------------------------------------------------------------------

dir.create(DOWNLOAD_DIR, showWarnings = FALSE)

download_month <- function(year, month) {
  out <- file.path(DOWNLOAD_DIR, sprintf("era5_%d_%02d.nc", year, month))
  if (file.exists(out)) return(out)

  wf_request(
    user    = wf_get_key(service = "cds"),
    request = list(
      dataset         = "derived-era5-single-levels-daily-statistics",
      product_type    = "reanalysis",
      variable        = "instantaneous_10m_wind_gust",
      year            = as.character(year),
      month           = sprintf("%02d", month),
      day             = sprintf("%02d", 1:31),
      daily_statistic = "daily_maximum",
      time_zone       = "utc+01:00",
      frequency       = "6_hourly",
      area            = AREA,
      format          = "netcdf"
    ),
    path   = DOWNLOAD_DIR,
    target = basename(out)
  )
  out
}

# Download
files <- mapply(download_month,
                rep(YEARS, each = 12),
                rep(1:12,  times = length(YEARS)))

# Merge along the time dimension using CDO if available, else ncdf4
if (nchar(Sys.which("cdo")) > 0) {
  system(paste("cdo mergetime", paste(files, collapse = " "), OUTPUT_FILE))
} else {
  # Manual merge with ncdf4
  datasets <- lapply(files, function(f) {
    nc  <- nc_open(f)
    on.exit(nc_close(nc))
    list(
      time = ncvar_get(nc, "time"),
      data = ncvar_get(nc, "i10fg"),
      lon  = ncvar_get(nc, "longitude"),
      lat  = ncvar_get(nc, "latitude"),
      tunits = ncatt_get(nc, "time", "units")$value
    )
  })

  lon  <- datasets[[1]]$lon
  lat  <- datasets[[1]]$lat
  time <- do.call(c, lapply(datasets, `[[`, "time"))
  data <- do.call(abind_time <- function(a, b) abind::abind(a, b, along = 3),
                  lapply(datasets, `[[`, "data"))

  d_lon  <- ncdim_def("longitude", "degrees_east",  lon)
  d_lat  <- ncdim_def("latitude",  "degrees_north", lat)
  d_time <- ncdim_def("time", datasets[[1]]$tunits,  time, unlim = TRUE)
  v_wind <- ncvar_def("i10fg", "m s-1", list(d_lon, d_lat, d_time),
                      missval = 1e20, longname = "Daily maximum wind gust")

  nc_out <- nc_create(OUTPUT_FILE, v_wind)
  ncvar_put(nc_out, v_wind, data)
  nc_close(nc_out)
}

message("Saved: ", OUTPUT_FILE)
message("Ready for: hm_read_netcdf('", OUTPUT_FILE, "', var = 'i10fg')")
