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
