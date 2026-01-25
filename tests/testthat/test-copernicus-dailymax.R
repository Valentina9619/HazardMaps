test_that("Copernicus dailymax NetCDF ingests with 1D lat/lon", {f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc", package = "HazardMaps")

                                                                 expect_true(nzchar(f))
                                                                 expect_true(file.exists(f))

                                                                 x <- hm_read_netcdf(f, var = "i10fg")
                                                                 x <- hm_standardize_coords(x)
                                                                 x <- hm_decode_time(x)
                                                                 X <- hm_to_matrix(x)

                                                                 expect_true(is.matrix(X))
                                                                 expect_equal(dim(X), c(10226L, 195L))

                                                                 coords <- attr(X, "coords")
                                                                 expect_true(is.data.frame(coords))
                                                                 expect_true(all(c("cell_id", "lat", "lon") %in% names(coords)))})
