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
