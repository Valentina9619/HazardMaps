# Extract grid geometry

Computes the bounding box, grid spacing, and cell count from the
coordinate table and stores them in `x$grid`.

## Usage

``` r
hm_get_grid_geometry(x)
```

## Arguments

- x:

  An `hm_hazard` object after \[hm_standardize_coords()\].

## Value

The same object with a new `$grid` element containing `nx`, `ny`, `dx`,
`dy`, `bbox`, `n_cells`, and `grid_type`.

## See also

\[hm_standardize_coords()\]
