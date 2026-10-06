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
