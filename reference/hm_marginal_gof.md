# Goodness-of-fit diagnostics for a single grid cell

Returns the KS statistic and QQ data for the best-fitting distribution
at a chosen cell. The KS test is a complementary diagnostic —
distribution selection is done by AIC/BIC in \[hm_fit_marginals()\].

## Usage

``` r
hm_marginal_gof(X, fits, cell_id)
```

## Arguments

- X:

  Numeric matrix passed to \[hm_fit_marginals()\].

- fits:

  data.frame from \[hm_fit_marginals()\].

- cell_id:

  Column index to diagnose.

## Value

A list with `cell_id`, `dist`, `dist_label`, `params`, `ks_statistic`,
`ks_pvalue`, `qq`, and `n_obs`.

## See also

\[hm_fit_marginals()\], \[hm_pit_transform()\]

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
gof  <- hm_marginal_gof(X, fits, cell_id = 1L)
gof$ks_statistic
} # }
```
