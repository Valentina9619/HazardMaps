test_that("Copernicus dailymax NetCDF ingests with 1D lat/lon", {f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc", package = "HazardMaps")

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
                                                                 expect_true(all(c("cell_id", "lat", "lon") %in% names(coords)))})


test_that("Copernicus dailymax spatial domain can be plotted", {f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc", package = "HazardMaps")

                                                                expect_true(nzchar(f))
                                                                expect_true(file.exists(f))

                                                                x <- hm_read_netcdf(f, var = "i10fg")
                                                                expect_equal(x$meta$grid_type, "regular")

                                                                x <- hm_standardize_coords(x)

                                                                p <- hm_plot_spatial_domain(x,
                                                                                            dataset = "copernicus",
                                                                                            bbox = c(lon_min = 10.75,
                                                                                                     lon_max = 14.00,
                                                                                                     lat_min = 43.75,
                                                                                                     lat_max = 47.25),
                                                                                            title = "Copernicus domain",
                                                                                            show_map = FALSE)

                                                                expect_true(inherits(p, "ggplot"))})
