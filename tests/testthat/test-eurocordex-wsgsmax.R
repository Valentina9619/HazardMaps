test_that("EuroCordex windgust dailymax NetCDF 1970 to 2005", {f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc", package = "HazardMaps")

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

                                                               if (all(c("rlat", "rlon") %in% names(coords))) {expect_true(is.numeric(coords$rlat))
                                                                                                               expect_true(is.numeric(coords$rlon))}})


test_that("EuroCordex spatial domain can be plotted", {f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc", package = "HazardMaps")

                                                       expect_true(nzchar(f))
                                                       expect_true(file.exists(f))

                                                       x <- hm_read_netcdf(f, var = "wsgsmax")
                                                       expect_equal(x$meta$grid_type, "rotated")

                                                       x <- hm_standardize_coords(x)

                                                       p <- hm_plot_spatial_domain(x,
                                                                                   dataset = "eurocordex",
                                                                                   bbox = c(lon_min = 6.5,
                                                                                            lon_max = 14.2,
                                                                                            lat_min = 43.7,
                                                                                            lat_max = 47.2),
                                                                                   title = "EURO-CORDEX domain",
                                                                                   show_map = FALSE)

                                                       expect_true(inherits(p, "ggplot"))})


test_that("EURO-CORDEX exceedance count can be plotted", {f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc", package = "HazardMaps")

                                                          expect_true(nzchar(f))
                                                          expect_true(file.exists(f))

                                                          x <- hm_read_netcdf(f, var = "wsgsmax")
                                                          x <- hm_standardize_coords(x)

                                                          p <- hm_plot_exceedance(x,
                                                                                  threshold = 25,
                                                                                  metric = "count",
                                                                                  dataset = "eurocordex",
                                                                                  bbox = c(lon_min = 6.5,
                                                                                           lon_max = 14.2,
                                                                                           lat_min = 43.7,
                                                                                           lat_max = 47.2),
                                                                                  title = "Wind-gust exceedance count",
                                                                                  show_map = FALSE)

                                                          expect_true(inherits(p, "ggplot"))})


test_that("EURO-CORDEX exceedance time plot works", {f <- system.file("extdata", "EuroCordex_wsgsmax_1970_2005.nc", package = "HazardMaps")

                                                     expect_true(nzchar(f))
                                                     expect_true(file.exists(f))

                                                     x <- hm_read_netcdf(f, var = "wsgsmax")
                                                     x <- hm_standardize_coords(x)

                                                     p_2y <- hm_plot_exceedance_time(x,
                                                                                     threshold = 25,
                                                                                     by = "2_years",
                                                                                     title = "EURO-CORDEX extreme-event days by 2-year period")

                                                     p_decade <- hm_plot_exceedance_time(x,
                                                                                         threshold = 25,
                                                                                         by = "decade",
                                                                                         title = "EURO-CORDEX extreme-event days by decade")

                                                     expect_true(inherits(p_2y, "ggplot"))
                                                     expect_true(inherits(p_decade, "ggplot"))})


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

  # EURO-CORDEX is a rotated grid that includes sea cells — those cells
  # have all-NA or all-zero wind-gust values and will have no successful
  # fit. We check that at least the majority of cells were fitted.
  n_fitted <- sum(!is.na(fits$best_dist))
  expect_gt(n_fitted, ncol(X) * 0.5)

  # Among fitted cells, Weibull must be the dominant distribution
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

  # Use the first successfully fitted cell for the diagnostic
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

  # U and Z must have the same dimensions as X
  expect_equal(dim(pit$U), dim(X))
  expect_equal(dim(pit$Z), dim(X))

  # Among transformed cells, U must be in (0, 1) and Z must be finite
  ok_cols <- which(!is.na(fits$best_dist))
  expect_true(all(pit$U[, ok_cols] >= 0 & pit$U[, ok_cols] <= 1,
                  na.rm = TRUE))
  expect_true(all(is.finite(pit$Z[, ok_cols]), na.rm = TRUE))

  # Skipped cells correspond to sea/no-data cells — their count must
  # match the number of cells with no successful fit
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
