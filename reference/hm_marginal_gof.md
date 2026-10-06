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
