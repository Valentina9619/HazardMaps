# Tests for the Copernicus daily maximum wind-gust dataset (regular grid).
# Covers the full pipeline: ingest → coordinates → time → matrix → event
# selection → event summary → visualization.

test_that("Copernicus dailymax NetCDF ingests with 1D lat/lon", {
  f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "i10fg")
  expect_equal(x$meta$grid_type, "regular")

  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  X <- hm_to_matrix(x)

  expect_true(is.matrix(X))
  expect_equal(dim(X), c(10226L, 195L))

  coords <- attr(X, "coords")
  expect_true(is.data.frame(coords))
  expect_true(all(c("cell_id", "lat", "lon") %in% names(coords)))
})


test_that("Copernicus dailymax spatial domain can be plotted", {
  f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "i10fg")
  expect_equal(x$meta$grid_type, "regular")

  x <- hm_standardize_coords(x)

  p <- hm_plot_spatial_domain(x,
                               dataset  = "copernicus",
                               bbox     = c(lon_min = 10.75, lon_max = 14.00,
                                            lat_min = 43.75, lat_max = 47.25),
                               title    = "Copernicus domain",
                               show_map = FALSE)

  expect_true(inherits(p, "ggplot"))
})


test_that("Copernicus exceedance count can be plotted", {
  f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "i10fg")
  x <- hm_standardize_coords(x)

  p <- hm_plot_exceedance(x,
                           threshold  = 25,
                           metric     = "count",
                           dataset    = "copernicus",
                           bbox       = c(lon_min = 11.5, lon_max = 15.0,
                                          lat_min = 44.0, lat_max = 47.0),
                           title      = "Wind-gust exceedance count",
                           show_map   = FALSE)

  expect_true(inherits(p, "ggplot"))
})


test_that("Copernicus exceedance time plot works", {
  f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "i10fg")
  x <- hm_standardize_coords(x)

  p_year <- hm_plot_exceedance_time(x,
                                     threshold = 25,
                                     by        = "year",
                                     title     = "Copernicus extreme-event days by year")

  p_5y <- hm_plot_exceedance_time(x,
                                   threshold = 25,
                                   by        = "5_years",
                                   title     = "Copernicus extreme-event days by 5-year period")

  expect_true(inherits(p_year, "ggplot"))
  expect_true(inherits(p_5y,   "ggplot"))
})


test_that("Copernicus extreme events can be selected by threshold", {
  f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "i10fg")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)

  x_ext <- hm_select_extreme_events(x, threshold = 25)

  expect_s3_class(x_ext, "hm_hazard")
  expect_lt(length(x_ext$time), length(x$time))
  expect_equal(x_ext$meta$filter$threshold,       25)
  expect_equal(x_ext$meta$filter$n_days_original, length(x$time))
  expect_equal(x_ext$meta$filter$n_days_retained, length(x_ext$time))

  # $data time dimension must match the number of retained dates
  d    <- dim(x_ext$data)
  tidx <- which(d == length(x_ext$time))
  expect_length(tidx, 1L)

  # Every retained day must have at least one cell >= threshold
  X       <- hm_to_matrix(x_ext)
  row_max <- apply(X, 1, max, na.rm = TRUE)
  expect_true(all(row_max >= 25))
})


test_that("Copernicus extreme-event days can be summarised", {
  f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "i10fg")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)

  x_ext <- hm_select_extreme_events(x, threshold = 25)
  out   <- hm_summarise_events(x_ext)

  expect_s3_class(out, "data.frame")
  expect_true(all(c("date", "n_cells_exceeded", "max_value") %in% names(out)))
  expect_equal(nrow(out), length(x_ext$time))
  expect_true(all(out$n_cells_exceeded >= 1L))
  expect_true(all(out$max_value >= 25, na.rm = TRUE))
})
