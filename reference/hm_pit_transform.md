# Apply the Probability Integral Transform

Transforms `X` to uniform margins using each cell's fitted CDF, then
optionally applies the probit transform to obtain centred Gaussian
scores.

## Usage

``` r
hm_pit_transform(X, fits, gaussian = TRUE)
```

## Arguments

- X:

  Numeric matrix from \[hm_fit_marginals()\].

- fits:

  A data frame returned by \[hm_fit_marginals()\], with one row per grid
  cell containing the best-fit distribution and parameters.

- gaussian:

  If `TRUE` (default), also returns centred Gaussian scores \\Z =
  \Phi^{-1}(U)\\.

## Value

A list with `U` (uniform matrix), `Z` (Gaussian scores or `NULL`),
`n_cells_ok`, and `n_cells_skipped`.

## See also

\[hm_fit_marginals()\], \[hm_marginal_gof()\]

## Examples

``` r
if (FALSE) { # \dontrun{
f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                 package = "HazardMaps")
x <- hm_read_netcdf(f, var = "i10fg")
x <- hm_standardize_coords(x)
x <- hm_decode_time(x)
x_ext <- hm_select_extreme_events(x, threshold = 25)
X    <- hm_to_matrix(x_ext)
fits <- hm_fit_marginals(X)
pit  <- hm_pit_transform(X, fits)
dim(pit$U)
} # }
```
