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

  d    <- dim(x_ext$data)
  tidx <- which(d == length(x_ext$time))
  expect_length(tidx, 1L)

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

  dist_counts <- table(fits$best_dist)
  expect_true("weibull" %in% names(dist_counts))
  expect_true(dist_counts["weibull"] == max(dist_counts))

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

  expect_true(is.list(pit))
  expect_true(all(c("U", "Z", "n_cells_ok", "n_cells_skipped") %in% names(pit)))

  expect_equal(dim(pit$U), dim(X))
  expect_equal(dim(pit$Z), dim(X))

  expect_true(all(pit$U >= 0 & pit$U <= 1, na.rm = TRUE))

  expect_true(all(is.finite(pit$Z), na.rm = TRUE))

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


test_that("Copernicus C-vine copula can be fitted and simulated", {
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

  U_small    <- pit$U[, 1:10]
  fits_small <- fits[1:10, ]

  cop <- hm_fit_copula(U_small, type = "cvine")

  expect_s3_class(cop, "hm_copula")
  expect_equal(cop$type,    "cvine")
  expect_equal(cop$n_cells, 10L)
  expect_equal(cop$n_obs,   nrow(U_small))

  scenarios <- hm_simulate_copula(cop, fits_small,
                                   n_scenarios = 50L, seed = 1L)

  expect_true(is.matrix(scenarios))
  expect_equal(dim(scenarios), c(50L, 10L))
  expect_true(all(scenarios >= 0, na.rm = TRUE))
  expect_equal(attr(scenarios, "type"), "cvine")
})


test_that("Copernicus D-vine copula can be fitted and simulated", {
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

  U_small    <- pit$U[, 1:10]
  fits_small <- fits[1:10, ]

  cop <- hm_fit_copula(U_small, type = "dvine")

  expect_s3_class(cop, "hm_copula")
  expect_equal(cop$type,    "dvine")
  expect_equal(cop$n_cells, 10L)

  scenarios <- hm_simulate_copula(cop, fits_small,
                                   n_scenarios = 50L, seed = 1L)

  expect_true(is.matrix(scenarios))
  expect_equal(dim(scenarios), c(50L, 10L))
  expect_true(all(scenarios >= 0, na.rm = TRUE))
  expect_equal(attr(scenarios, "type"), "dvine")
})


test_that("Copernicus Gaussian copula can be fitted and simulated", {
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

  U_small    <- pit$U[, 1:10]
  fits_small <- fits[1:10, ]

  cop <- hm_fit_copula(U_small, type = "gaussian")

  expect_s3_class(cop, "hm_copula")
  expect_equal(cop$type,    "gaussian")
  expect_equal(cop$n_cells, 10L)

  R <- cop$model$R
  expect_true(isSymmetric(R))
  expect_true(all(eigen(R)$values > 0))

  scenarios <- hm_simulate_copula(cop, fits_small,
                                   n_scenarios = 50L, seed = 1L)

  expect_true(is.matrix(scenarios))
  expect_equal(dim(scenarios), c(50L, 10L))
  expect_true(all(scenarios >= 0, na.rm = TRUE))
  expect_equal(attr(scenarios, "type"), "gaussian")
})


test_that("Copernicus spatial correlation model can be fitted", {
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

  model <- hm_fit_correlation_model(corg)

  expect_s3_class(model, "hm_corrmodel")
  expect_true(!is.null(model$best_family))
  expect_true(model$best_family %in%
                c("matern", "exponential", "gaussian", "spherical"))
  expect_true(all(is.finite(model$best_params)))
  expect_true(model$best_rmse < 0.5)

  expect_s3_class(model$all_fits, "data.frame")
  expect_equal(nrow(model$all_fits), 4L)
  expect_true(all(model$all_fits$converged))
})


test_that("Copernicus KLE scenarios have correct structure and range", {
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
  model  <- hm_fit_correlation_model(corg)

  scenarios <- hm_simulate_sf_kle(model, D, fits[1:20, ],
                                   n_scenarios = 50L,
                                   var_retain  = 0.95,
                                   seed        = 42L)

  expect_true(is.matrix(scenarios))
  expect_equal(dim(scenarios), c(50L, 20L))
  expect_true(all(scenarios > 0, na.rm = TRUE))

  expect_true(!is.null(attr(scenarios, "k_modes")))
  expect_true(!is.null(attr(scenarios, "var_explained")))
  expect_lte(attr(scenarios, "k_modes"), 20L)
  expect_gte(attr(scenarios, "var_explained"), 0.95)
})
