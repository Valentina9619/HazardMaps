# Summarise retained extreme-event days

Returns a data.frame with one row per time step: date, number of cells
that met the threshold, and domain-wide maximum value.

## Usage

``` r
hm_summarise_events(x, threshold = NULL)
```

## Arguments

- x:

  An `hm_hazard` object after \[hm_standardize_coords()\] and
  \[hm_decode_time()\].

- threshold:

  Threshold used to count exceeding cells. If `NULL`, read from
  `x$meta$filter$threshold`.

## Value

A data.frame with columns `date`, `n_cells_exceeded`, and `max_value`.

## See also

\[hm_select_extreme_events()\]

## Examples

``` r
f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                 package = "HazardMaps")
x <- hm_read_netcdf(f, var = "i10fg")
x <- hm_standardize_coords(x)
x <- hm_decode_time(x)
x_ext <- hm_select_extreme_events(x, threshold = 25)
ev <- hm_summarise_events(x_ext)
head(ev)
#>         date n_cells_exceeded max_value
#> 1 1994-01-01               24  30.32205
#> 2 1994-01-06                2  25.50521
#> 3 1994-01-07                8  27.78774
#> 4 1994-01-20                1  25.16436
#> 5 1994-01-21                1  25.72410
#> 6 1994-01-26               10  30.36147
```
