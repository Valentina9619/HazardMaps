# Select extreme-event days by threshold

Retains only the time steps in which at least one grid cell reaches or
exceeds `threshold`. Call after \[hm_decode_time()\] and before
\[hm_to_matrix()\].

## Usage

``` r
hm_select_extreme_events(x, threshold)
```

## Arguments

- x:

  An `hm_hazard` object after \[hm_standardize_coords()\] and
  \[hm_decode_time()\].

- threshold:

  A single finite numeric value.

## Value

The same object with `$data` and `$time` restricted to extreme-event
days. `$meta$filter` records the threshold and day counts for
reproducibility.

## See also

\[hm_summarise_events()\]
