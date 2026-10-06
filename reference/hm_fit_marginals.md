# Fit marginal distributions to each grid cell

Fits candidate parametric distributions to each column of `X` by maximum
likelihood and selects the best fit using AIC or BIC, following the
methodology in Clavijo Mesa et al. (2026).

## Usage

``` r
hm_fit_marginals(X, dists = .hm_supported_dists, criterion = "bic")
```

## Arguments

- X:

  A numeric matrix (\\T \times N\\) from \[hm_to_matrix()\]. All values
  must be strictly positive.

- dists:

  Candidate families. Any subset of
  `c("weibull", "gamma", "lnorm", "gumbel", "tnorm", "tlaplace")`.

- criterion:

  Selection criterion: `"bic"` (default) or `"aic"`.

## Value

A data.frame with one row per cell and columns `cell_id`, `best_dist`,
`best_dist_label`, `params`, `loglik`, `aic`, `bic`, `n_obs`,
`n_failed`. The attribute `"fits_full"` stores raw per-cell fit results
for use in \[hm_marginal_gof()\].

## See also

\[hm_marginal_gof()\], \[hm_pit_transform()\]
