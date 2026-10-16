# Plot extreme-event days over time

Bar chart of extreme-event days aggregated by year, two-year, five-year
or decade periods. A day is counted once when at least one cell exceeds
`threshold`.

## Usage

``` r
hm_plot_exceedance_time(
  x,
  threshold,
  by = "year",
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  bar_fill = "grey70",
  bar_colour = "grey30"
)
```

## Arguments

- x:

  An `hm_hazard` object after \[hm_standardize_coords()\].

- threshold:

  Numeric threshold defining exceedances.

- by:

  Aggregation period: `"year"`, `"2_years"`, `"5_years"`, or `"decade"`.

- title:

  Optional plot title.

- xlab:

  Optional x-axis label.

- ylab:

  Optional y-axis label.

- bar_fill:

  Bar fill colour.

- bar_colour:

  Bar border colour.

## Value

A `ggplot` object.

## Examples

``` r
f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                 package = "HazardMaps")
x <- hm_read_netcdf(f, var = "i10fg")
x <- hm_standardize_coords(x)
x <- hm_decode_time(x)
hm_plot_exceedance_time(x, threshold = 25, by = "year")
```
