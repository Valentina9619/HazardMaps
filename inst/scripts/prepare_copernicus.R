# Download ERA5 daily-maximum wind-gust files from Copernicus CDS and merge
# them into a single NetCDF ready for hm_read_netcdf().
#
# Requirements: ecmwfr, ncdf4
#   install.packages(c("ecmwfr", "ncdf4"))
#
# Before running: register at https://cds.climate.copernicus.eu and store
#   your credentials once with ecmwfr::wf_set_key()

library(ecmwfr)
library(ncdf4)

# Configuration -----------------------------------------------------------
YEARS        <- 1994:2021
AREA         <- c(47.25, 10.75, 43.75, 14.25)  # N / W / S / E bounding box
DOWNLOAD_DIR <- "era5_monthly"
OUTPUT_FILE  <- "Copernicus_dailymax_1994_2021.nc"
# --------------------------------------------------------------------------

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


# Download one file per month (skips files already on disk)
files <- mapply(download_month,
                rep(YEARS, each = 12),
                rep(1:12,  times = length(YEARS)))

# Merge all monthly files into one NetCDF along the time dimension
message("Merging ", length(files), " monthly files ...")

nc1    <- nc_open(files[[1]])
lon    <- ncvar_get(nc1, "longitude")
lat    <- ncvar_get(nc1, "latitude")
tunits <- ncatt_get(nc1, "time", "units")$value
nc_close(nc1)

all_time <- unlist(lapply(files, function(f) {
  nc <- nc_open(f); on.exit(nc_close(nc))
  ncvar_get(nc, "time")
}))

d_lon  <- ncdim_def("longitude", "degrees_east",  lon)
d_lat  <- ncdim_def("latitude",  "degrees_north", lat)
d_time <- ncdim_def("time", tunits, all_time, unlim = TRUE)
v_wind <- ncvar_def("i10fg", "m s-1", list(d_lon, d_lat, d_time),
                    missval = 1e20, longname = "Daily maximum wind gust")

nc_out <- nc_create(OUTPUT_FILE, v_wind)

t_start <- 1L
for (f in files) {
  nc   <- nc_open(f)
  dat  <- ncvar_get(nc, "i10fg")
  nt   <- dim(dat)[3]
  ncvar_put(nc_out, v_wind, dat,
            start = c(1L, 1L, t_start),
            count = c(-1L, -1L, nt))
  t_start <- t_start + nt
  nc_close(nc)
}

nc_close(nc_out)
message("Saved: ", OUTPUT_FILE)
message("Ready for: hm_read_netcdf('", OUTPUT_FILE, "', var = 'i10fg')")
