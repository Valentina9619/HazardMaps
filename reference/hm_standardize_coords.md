# Standardize spatial coordinates

Converts raw NetCDF coordinates into a canonical data.frame with one row
per grid cell.

## Usage

``` r
hm_standardize_coords(x)
```

## Arguments

- x:

  An `hm_hazard` object returned by \[hm_read_netcdf()\].

## Value

The same object with `x$coords` replaced by a data.frame containing
`cell_id`, `lat`, and `lon`.

## Details

Regular grids supply `lat` and `lon` as 1D vectors; all combinations are
expanded with
[`expand.grid()`](https://rdrr.io/r/base/expand.grid.html). Rotated or
curvilinear grids supply `lat` and `lon` as 2D matrices, with rotated
indices appended when available.

## See also

\[hm_read_netcdf()\], \[hm_get_grid_geometry()\]

## Examples

``` r
f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                 package = "HazardMaps")
x <- hm_read_netcdf(f, var = "i10fg")
x <- hm_standardize_coords(x)
head(x$coords)
#>   cell_id   lat  lon
#> 1       1 47.00 11.5
#> 2       2 46.75 11.5
#> 3       3 46.50 11.5
#> 4       4 46.25 11.5
#> 5       5 46.00 11.5
#> 6       6 45.75 11.5
```
