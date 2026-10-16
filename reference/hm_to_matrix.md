# Convert hazard array to a time-by-cell matrix

Reshapes `x$data` into a \\T \times N\\ numeric matrix where rows are
time steps and columns are grid cells.

## Usage

``` r
hm_to_matrix(x)
```

## Arguments

- x:

  An `hm_hazard` object after \[hm_standardize_coords()\] and
  \[hm_decode_time()\].

## Value

A numeric matrix with attributes `time` and `coords`.

## See also

\[hm_decode_time()\], \[hm_standardize_coords()\]

## Examples

``` r
f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                 package = "HazardMaps")
x <- hm_read_netcdf(f, var = "i10fg")
x <- hm_standardize_coords(x)
x <- hm_decode_time(x)
X <- hm_to_matrix(x)
dim(X)
#> [1] 10226   195
```
