# Remove sea and no-data cells from a hazard matrix

Drops columns where the fraction of finite positive values falls below
`min_valid`. Needed for rotated-grid datasets such as EURO-CORDEX where
sea cells produce all-NA columns.

## Usage

``` r
hm_filter_land_cells(X, min_valid = 0.8)
```

## Arguments

- X:

  A numeric matrix as returned by \[hm_to_matrix()\].

- min_valid:

  Minimum fraction of valid observations required to keep a column.
  Defaults to `0.8`.

## Value

A matrix with only land columns retained. Attributes `land_idx`,
`coords`, and `time` are attached.

## See also

\[hm_to_matrix()\], \[hm_fit_marginals()\]
