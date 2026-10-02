# Tests for the EURO-CORDEX daily maximum wind-gust dataset (rotated grid).
# Covers the full pipeline: ingest → coordinates → time → matrix → event
# selection → event summary → visualization.

test_that("EuroCordex windgust dailymax NetCDF 1970 to 2005", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
  expect_equal(x$meta$grid_type, "rotated")

  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  X <- hm_to_matrix(x)

  expect_true(is.matrix(X))
  expect_equal(dim(X), c(13149L, 1872L))

  coords <- attr(X, "coords")
  expect_true(is.data.frame(coords))
  expect_true(all(c("cell_id", "lat", "lon") %in% names(coords)))

  if (all(c("rlat", "rlon") %in% names(coords))) {
    expect_true(is.numeric(coords$rlat))
    expect_true(is.numeric(coords$rlon))
  }
})


test_that("EuroCordex spatial domain can be plotted", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
  expect_equal(x$meta$grid_type, "rotated")

  x <- hm_standardize_coords(x)

  p <- hm_plot_spatial_domain(x,
                               dataset  = "eurocordex",
                               bbox     = c(lon_min = 6.5,  lon_max = 14.2,
                                            lat_min = 43.7, lat_max = 47.2),
                               title    = "EURO-CORDEX domain",
                               show_map = FALSE)

  expect_true(inherits(p, "ggplot"))
})


test_that("EURO-CORDEX exceedance count can be plotted", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
  x <- hm_standardize_coords(x)

  p <- hm_plot_exceedance(x,
                           threshold  = 25,
                           metric     = "count",
                           dataset    = "eurocordex",
                           bbox       = c(lon_min = 6.5,  lon_max = 14.2,
                                          lat_min = 43.7, lat_max = 47.2),
                           title      = "Wind-gust exceedance count",
                           show_map   = FALSE)

  expect_true(inherits(p, "ggplot"))
})


test_that("EURO-CORDEX exceedance time plot works", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
  x <- hm_standardize_coords(x)

  p_2y <- hm_plot_exceedance_time(x,
                                   threshold = 25,
                                   by        = "2_years",
                                   title     = "EURO-CORDEX extreme-event days by 2-year period")

  p_decade <- hm_plot_exceedance_time(x,
                                       threshold = 25,
                                       by        = "decade",
                                       title     = "EURO-CORDEX extreme-event days by decade")

  expect_true(inherits(p_2y,    "ggplot"))
  expect_true(inherits(p_decade, "ggplot"))
})


test_that("EURO-CORDEX extreme events can be selected by threshold", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
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


test_that("EURO-CORDEX extreme-event days can be summarised", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
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
