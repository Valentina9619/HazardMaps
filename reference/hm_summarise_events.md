# Summarise retained extreme-event days

Returns a data.frame with one row per time step: date, number of cells
that met the threshold, and domain-wide maximum value.

## Usage

``` r
hm_summarise_events(x, threshold = NULL)
```

## Arguments

- x:

  An `hm_hazard` object, typically from \[hm_select_extreme_events()\].

- threshold:

  Threshold used to count exceeding cells. If `NULL`, read from
  `x$meta$filter$threshold`.

## Value

A data.frame with columns `date`, `n_cells_exceeded`, and `max_value`.

## See also

\[hm_select_extreme_events()\]
