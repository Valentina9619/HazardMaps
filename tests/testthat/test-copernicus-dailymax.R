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


test_that("Copernicus marginal distributions can be fitted", {
  f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "i10fg")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x <- hm_select_extreme_events(x, threshold = 25)
  X <- hm_to_matrix(x)

  fits <- hm_fit_marginals(X, criterion = "bic")

  expect_s3_class(fits, "data.frame")
  expect_equal(nrow(fits), ncol(X))
  expect_true(all(c("cell_id", "best_dist", "aic", "bic",
                    "params", "n_obs", "n_failed") %in% names(fits)))

  # Weibull must be the dominant distribution for wind-gust data
  dist_counts <- table(fits$best_dist)
  expect_true("weibull" %in% names(dist_counts))
  expect_true(dist_counts["weibull"] == max(dist_counts))

  # Every cell must have enough observations and a successful fit
  expect_true(all(fits$n_obs > 0L))
  expect_true(all(!is.na(fits$best_dist)))
})


test_that("Copernicus goodness-of-fit diagnostic works", {
  f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "i10fg")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x <- hm_select_extreme_events(x, threshold = 25)
  X <- hm_to_matrix(x)

  fits <- hm_fit_marginals(X)
  gof  <- hm_marginal_gof(X, fits, cell_id = 1L)

  expect_true(is.list(gof))
  expect_true(all(c("cell_id", "dist", "params",
                    "ks_statistic", "ks_pvalue", "qq",
                    "n_obs") %in% names(gof)))

  expect_equal(gof$cell_id, 1L)
  expect_true(is.numeric(gof$ks_statistic) && gof$ks_statistic >= 0)
  expect_true(is.data.frame(gof$qq))
  expect_true(all(c("theoretical", "empirical") %in% names(gof$qq)))

  # A good fit: KS statistic should be small (< 0.2 for wind-gust data)
  expect_lt(gof$ks_statistic, 0.2)
})


test_that("Copernicus PIT transform produces valid uniform and Gaussian matrices", {
  f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "i10fg")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x <- hm_select_extreme_events(x, threshold = 25)
  X <- hm_to_matrix(x)

  fits <- hm_fit_marginals(X)
  pit  <- hm_pit_transform(X, fits)

  # Output structure
  expect_true(is.list(pit))
  expect_true(all(c("U", "Z", "n_cells_ok", "n_cells_skipped") %in% names(pit)))

  # U and Z must have the same dimensions as X
  expect_equal(dim(pit$U), dim(X))
  expect_equal(dim(pit$Z), dim(X))

  # U must be in (0, 1) — it is the uniform PIT output
  expect_true(all(pit$U >= 0 & pit$U <= 1, na.rm = TRUE))

  # Z can be negative — it is the Gaussian anamorphosis of U
  expect_true(all(is.finite(pit$Z), na.rm = TRUE))

  # All cells must have been transformed successfully
  expect_equal(pit$n_cells_ok,      ncol(X))
  expect_equal(pit$n_cells_skipped, 0L)
})


test_that("Copernicus best-fit distribution map can be plotted", {
  f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "i10fg")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x_ext <- hm_select_extreme_events(x, threshold = 25)
  X     <- hm_to_matrix(x_ext)
  fits  <- hm_fit_marginals(X)

  p <- hm_plot_best_dist(x_ext, fits, show_map = FALSE)

  expect_true(inherits(p, "ggplot"))
})


test_that("Copernicus marginal fit diagnostic can be plotted", {
  f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "i10fg")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x_ext <- hm_select_extreme_events(x, threshold = 25)
  X     <- hm_to_matrix(x_ext)
  fits  <- hm_fit_marginals(X)

  result <- hm_plot_marginal_fit(X, fits, cell_id = 1L, x = x_ext)

  # Returns a patchwork object or a list of ggplots depending on
  # whether patchwork is installed
  expect_true(inherits(result, "ggplot") || is.list(result))
})


test_that("Copernicus distribution frequency chart can be plotted", {
  f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "i10fg")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x_ext <- hm_select_extreme_events(x, threshold = 25)
  X     <- hm_to_matrix(x_ext)
  fits  <- hm_fit_marginals(X)

  p <- hm_plot_dist_frequency(fits)

  expect_true(inherits(p, "ggplot"))
})


test_that("Copernicus land cell filter retains all cells (regular grid)", {
  f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "i10fg")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x <- hm_select_extreme_events(x, threshold = 25)
  X <- hm_to_matrix(x)

  X_land <- hm_filter_land_cells(X)

  # Copernicus is a land-only regular grid — all cells must be retained
  expect_equal(ncol(X_land), ncol(X))
  expect_true(!is.null(attr(X_land, "land_idx")))
  expect_equal(length(attr(X_land, "land_idx")), ncol(X))
  expect_equal(nrow(X_land), nrow(X))
})


test_that("Copernicus distance matrix has correct structure", {
  f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "i10fg")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x <- hm_select_extreme_events(x, threshold = 25)
  X <- hm_to_matrix(x)

  X_land <- hm_filter_land_cells(X)
  coords <- attr(X_land, "coords")

  # Use a small subset for speed
  D <- hm_distance_matrix(coords[1:20, ])

  expect_true(is.matrix(D))
  expect_equal(dim(D), c(20L, 20L))
  expect_true(all(diag(D) == 0))
  expect_true(all(D >= 0))
  expect_true(isSymmetric(D))
  expect_gt(max(D), 0)
})


test_that("Copernicus empirical Kendall matrix has correct structure", {
  f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "i10fg")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x <- hm_select_extreme_events(x, threshold = 25)
  X <- hm_to_matrix(x)

  X_land <- hm_filter_land_cells(X)
  fits   <- hm_fit_marginals(X_land)
  pit    <- hm_pit_transform(X_land, fits)
  tau    <- hm_empirical_kendall(pit$Z, cell_idx = 1:20)

  expect_true(is.matrix(tau))
  expect_equal(dim(tau), c(20L, 20L))
  expect_true(all(diag(tau) == 1))
  expect_true(isSymmetric(tau))
  expect_true(all(tau >= -1 & tau <= 1, na.rm = TRUE))
})


test_that("Copernicus empirical correlogram decreases with distance", {
  f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "i10fg")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x <- hm_select_extreme_events(x, threshold = 25)
  X <- hm_to_matrix(x)

  X_land <- hm_filter_land_cells(X)
  fits   <- hm_fit_marginals(X_land)
  pit    <- hm_pit_transform(X_land, fits)
  tau    <- hm_empirical_kendall(pit$Z, cell_idx = 1:20)
  D      <- hm_distance_matrix(attr(X_land, "coords")[1:20, ])
  corg   <- hm_empirical_correlogram(tau, D, n_bins = 10)

  expect_s3_class(corg, "data.frame")
  expect_true(all(c("bin_center", "mean_tau", "n_pairs",
                    "bin_lo", "bin_hi") %in% names(corg)))
  expect_true(all(corg$n_pairs > 0L))
  expect_true(all(corg$bin_center > 0))

  # Correlation at the shortest distances must be higher than at long range
  expect_gt(corg$mean_tau[1], corg$mean_tau[nrow(corg)])
})


test_that("Copernicus cluster size output is consistent", {
  f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "i10fg")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x <- hm_select_extreme_events(x, threshold = 25)
  X <- hm_to_matrix(x)

  X_land <- hm_filter_land_cells(X)
  cs     <- hm_cluster_size(X_land, threshold = 25,
                             time = attr(X_land, "time"))

  expect_s3_class(cs, "data.frame")
  expect_true(all(c("date", "cluster_size",
                    "cluster_fraction") %in% names(cs)))
  expect_equal(nrow(cs), nrow(X_land))

  # Every retained day must have at least one cell >= threshold
  # (because X_land comes from hm_select_extreme_events)
  expect_true(all(cs$cluster_size >= 1L))
  expect_true(all(cs$cluster_fraction > 0 & cs$cluster_fraction <= 1))
})


test_that("Copernicus correlogram can be plotted", {
  f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "i10fg")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x <- hm_select_extreme_events(x, threshold = 25)
  X <- hm_to_matrix(x)

  X_land <- hm_filter_land_cells(X)
  fits   <- hm_fit_marginals(X_land)
  pit    <- hm_pit_transform(X_land, fits)
  tau    <- hm_empirical_kendall(pit$Z, cell_idx = 1:20)
  D      <- hm_distance_matrix(attr(X_land, "coords")[1:20, ])
  corg   <- hm_empirical_correlogram(tau, D, n_bins = 10)

  p <- hm_plot_correlogram(corg,
         title = "Empirical correlogram — Copernicus wind gust")

  expect_true(inherits(p, "ggplot"))
})
