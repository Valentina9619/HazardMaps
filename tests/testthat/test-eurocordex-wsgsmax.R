test_that("EURO-CORDEX wsgsmax NetCDF ingests with 2D lat/lon", {
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


test_that("EURO-CORDEX spatial domain can be plotted", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
  expect_equal(x$meta$grid_type, "rotated")

  x <- hm_standardize_coords(x)

  p <- hm_plot_spatial_domain(x,
                               dataset  = "eurocordex",
                               bbox     = c(lon_min = 6.5, lon_max = 14.2,
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
                           threshold = 25,
                           metric    = "count",
                           dataset   = "eurocordex",
                           bbox      = c(lon_min = 6.5, lon_max = 14.2,
                                         lat_min = 43.7, lat_max = 47.2),
                           title     = "Wind-gust exceedance count",
                           show_map  = FALSE)

  expect_true(inherits(p, "ggplot"))
})


test_that("EURO-CORDEX exceedance time plot works", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
  x <- hm_standardize_coords(x)

  p_2y <- hm_plot_exceedance_time(x, threshold = 25, by = "2_years",
                                   title = "EURO-CORDEX extreme-event days by 2-year period")

  p_decade <- hm_plot_exceedance_time(x, threshold = 25, by = "decade",
                                       title = "EURO-CORDEX extreme-event days by decade")

  expect_true(inherits(p_2y,     "ggplot"))
  expect_true(inherits(p_decade, "ggplot"))
})


test_that("EURO-CORDEX marginal distributions can be fitted", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x <- hm_select_extreme_events(x, threshold = 25)
  X <- hm_to_matrix(x)

  fits <- hm_fit_marginals(X, criterion = "bic")

  expect_s3_class(fits, "data.frame")
  expect_equal(nrow(fits), ncol(X))
  expect_true(all(c("cell_id", "best_dist", "aic", "bic",
                    "params", "n_obs", "n_failed") %in% names(fits)))

  n_fitted <- sum(!is.na(fits$best_dist))
  expect_gt(n_fitted, ncol(X) * 0.5)

  dist_counts <- table(fits$best_dist[!is.na(fits$best_dist)])
  expect_true("weibull" %in% names(dist_counts))
  expect_true(dist_counts["weibull"] == max(dist_counts))
})


test_that("EURO-CORDEX goodness-of-fit diagnostic works", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x <- hm_select_extreme_events(x, threshold = 25)
  X <- hm_to_matrix(x)

  fits <- hm_fit_marginals(X)

  fitted_cell <- which(!is.na(fits$best_dist))[1L]
  gof <- hm_marginal_gof(X, fits, cell_id = fitted_cell)

  expect_true(is.list(gof))
  expect_true(all(c("cell_id", "dist", "params",
                    "ks_statistic", "ks_pvalue", "qq",
                    "n_obs") %in% names(gof)))

  expect_equal(gof$cell_id, fitted_cell)
  expect_true(is.numeric(gof$ks_statistic) && gof$ks_statistic >= 0)
  expect_true(is.data.frame(gof$qq))
  expect_true(all(c("theoretical", "empirical") %in% names(gof$qq)))
  expect_lt(gof$ks_statistic, 0.2)
})


test_that("EURO-CORDEX PIT transform produces valid uniform and Gaussian matrices", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
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

  ok_cols <- which(!is.na(fits$best_dist))
  expect_true(all(pit$U[, ok_cols] >= 0 & pit$U[, ok_cols] <= 1,
                  na.rm = TRUE))
  expect_true(all(is.finite(pit$Z[, ok_cols]), na.rm = TRUE))

  expect_equal(pit$n_cells_skipped, sum(is.na(fits$best_dist)))
  expect_equal(pit$n_cells_ok + pit$n_cells_skipped, ncol(X))
})


test_that("EURO-CORDEX best-fit distribution map can be plotted", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x_ext <- hm_select_extreme_events(x, threshold = 25)
  X     <- hm_to_matrix(x_ext)
  fits  <- hm_fit_marginals(X)

  p <- hm_plot_best_dist(x_ext, fits, show_map = FALSE)

  expect_true(inherits(p, "ggplot"))
})


test_that("EURO-CORDEX marginal fit diagnostic can be plotted", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x_ext <- hm_select_extreme_events(x, threshold = 25)
  X     <- hm_to_matrix(x_ext)
  fits  <- hm_fit_marginals(X)

  # Use the first successfully fitted cell
  fitted_cell <- which(!is.na(fits$best_dist))[1L]
  result <- hm_plot_marginal_fit(X, fits, cell_id = fitted_cell, x = x_ext)

  expect_true(inherits(result, "ggplot") || is.list(result))
})


test_that("EURO-CORDEX distribution frequency chart can be plotted", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x_ext <- hm_select_extreme_events(x, threshold = 25)
  X     <- hm_to_matrix(x_ext)
  fits  <- hm_fit_marginals(X)

  p <- hm_plot_dist_frequency(fits)

  expect_true(inherits(p, "ggplot"))
})


test_that("EURO-CORDEX land cell filter removes sea cells", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x <- hm_select_extreme_events(x, threshold = 25)
  X <- hm_to_matrix(x)

  X_land <- hm_filter_land_cells(X)

  # EURO-CORDEX has sea cells — land filter must remove some
  expect_lt(ncol(X_land), ncol(X))
  expect_true(!is.null(attr(X_land, "land_idx")))
  expect_equal(length(attr(X_land, "land_idx")), ncol(X_land))
  expect_equal(nrow(X_land), nrow(X))

  # coords attribute must be restricted to land cells
  coords_land <- attr(X_land, "coords")
  expect_equal(nrow(coords_land), ncol(X_land))
})


test_that("EURO-CORDEX distance matrix has correct structure", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x <- hm_select_extreme_events(x, threshold = 25)
  X <- hm_to_matrix(x)

  X_land <- hm_filter_land_cells(X)
  coords <- attr(X_land, "coords")

  D <- hm_distance_matrix(coords[1:20, ])

  expect_true(is.matrix(D))
  expect_equal(dim(D), c(20L, 20L))
  expect_true(all(diag(D) == 0))
  expect_true(all(D >= 0))
  expect_true(isSymmetric(D))
  expect_gt(max(D), 0)
})


test_that("EURO-CORDEX empirical Kendall matrix has correct structure", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
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


test_that("EURO-CORDEX empirical correlogram has correct structure", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
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
  expect_gt(corg$mean_tau[1], corg$mean_tau[nrow(corg)])
})


test_that("EURO-CORDEX cluster size output is consistent", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
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


test_that("EURO-CORDEX correlogram can be plotted", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
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
         title = "Empirical correlogram — EURO-CORDEX wind gust")

  expect_true(inherits(p, "ggplot"))
})


test_that("EURO-CORDEX C-vine copula can be fitted and simulated", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x <- hm_select_extreme_events(x, threshold = 25)
  X <- hm_to_matrix(x)

  X_land <- hm_filter_land_cells(X)
  fits   <- hm_fit_marginals(X_land)
  pit    <- hm_pit_transform(X_land, fits)

  fitted_idx <- which(!is.na(fits$best_dist))[1:10]
  U_small    <- pit$U[, fitted_idx]
  fits_small <- fits[fitted_idx, ]

  cop <- hm_fit_copula(U_small, type = "cvine")

  expect_s3_class(cop, "hm_copula")
  expect_equal(cop$type,    "cvine")
  expect_equal(cop$n_cells, 10L)

  scenarios <- hm_simulate_copula(cop, fits_small,
                                   n_scenarios = 50L, seed = 1L)

  expect_true(is.matrix(scenarios))
  expect_equal(dim(scenarios), c(50L, 10L))
  expect_true(all(scenarios >= 0, na.rm = TRUE))
  expect_equal(attr(scenarios, "type"), "cvine")
})


test_that("EURO-CORDEX Gaussian copula can be fitted and simulated", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x <- hm_select_extreme_events(x, threshold = 25)
  X <- hm_to_matrix(x)

  X_land <- hm_filter_land_cells(X)
  fits   <- hm_fit_marginals(X_land)
  pit    <- hm_pit_transform(X_land, fits)

  fitted_idx <- which(!is.na(fits$best_dist))[1:10]
  U_small    <- pit$U[, fitted_idx]
  fits_small <- fits[fitted_idx, ]

  cop <- hm_fit_copula(U_small, type = "gaussian")

  expect_s3_class(cop, "hm_copula")
  expect_equal(cop$type,    "gaussian")

  scenarios <- hm_simulate_copula(cop, fits_small,
                                   n_scenarios = 50L, seed = 1L)

  expect_true(is.matrix(scenarios))
  expect_equal(dim(scenarios), c(50L, 10L))
  expect_true(all(scenarios >= 0, na.rm = TRUE))
})


test_that("EURO-CORDEX spatial correlation model can be fitted", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x <- hm_select_extreme_events(x, threshold = 25)
  X <- hm_to_matrix(x)

  X_land     <- hm_filter_land_cells(X)
  fits       <- hm_fit_marginals(X_land)
  pit        <- hm_pit_transform(X_land, fits)
  fitted_idx <- which(!is.na(fits$best_dist))[1:20]
  tau        <- hm_empirical_kendall(pit$Z, cell_idx = fitted_idx)
  D          <- hm_distance_matrix(attr(X_land, "coords")[fitted_idx, ])
  corg       <- hm_empirical_correlogram(tau, D, n_bins = 10)

  model <- hm_fit_correlation_model(corg)

  expect_s3_class(model, "hm_corrmodel")
  expect_true(model$best_family %in%
                c("matern", "exponential", "gaussian", "spherical"))
  expect_true(all(is.finite(model$best_params)))
  expect_true(model$best_rmse < 0.5)
  expect_equal(nrow(model$all_fits), 4L)
})


test_that("EURO-CORDEX KLE scenarios have correct structure and range", {
  f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc",
                   package = "HazardMaps")

  expect_true(nzchar(f))
  expect_true(file.exists(f))

  x <- hm_read_netcdf(f, var = "wsgsmax")
  x <- hm_standardize_coords(x)
  x <- hm_decode_time(x)
  x <- hm_select_extreme_events(x, threshold = 25)
  X <- hm_to_matrix(x)

  X_land     <- hm_filter_land_cells(X)
  fits       <- hm_fit_marginals(X_land)
  pit        <- hm_pit_transform(X_land, fits)
  fitted_idx <- which(!is.na(fits$best_dist))[1:20]
  tau        <- hm_empirical_kendall(pit$Z, cell_idx = fitted_idx)
  D          <- hm_distance_matrix(attr(X_land, "coords")[fitted_idx, ])
  corg       <- hm_empirical_correlogram(tau, D, n_bins = 10)
  model      <- hm_fit_correlation_model(corg)

  scenarios <- hm_simulate_sf_kle(model, D, fits[fitted_idx, ],
                                   n_scenarios = 50L,
                                   var_retain  = 0.95,
                                   seed        = 42L)

  expect_true(is.matrix(scenarios))
  expect_equal(dim(scenarios), c(50L, 20L))
  expect_true(all(scenarios > 0, na.rm = TRUE))
  expect_lte(attr(scenarios, "k_modes"), 20L)
  expect_gte(attr(scenarios, "var_explained"), 0.95)
})
