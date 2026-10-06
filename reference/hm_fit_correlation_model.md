# Fit a parametric spatial correlation model to the empirical correlogram

Fits Matérn, Exponential, Gaussian, and/or Spherical covariance families
to the correlogram from \[hm_empirical_correlogram()\] by weighted least
squares (Eq. 6–7 of Clavijo Mesa et al., 2026, SF paper). Selects the
best family by RMSE, matching Table 1 of that paper.

## Usage

``` r
hm_fit_correlation_model(
  corg,
  families = c("matern", "exponential", "gaussian", "spherical"),
  max_dist = NULL
)
```

## Arguments

- corg:

  data.frame from \[hm_empirical_correlogram()\].

- families:

  Covariance families to fit. Any subset of
  `c("matern", "exponential", "gaussian", "spherical")`. Defaults to all
  four.

- max_dist:

  Maximum distance in km to include. Defaults to the full correlogram
  range.

## Value

An `hm_corrmodel` list with `best_family`, `best_params`, `best_rmse`,
`all_fits`, and `corg`.

## References

Clavijo Mesa, M.V., Broggi, M., Di Maio, F. & Zio, E. (2026).
*International Journal of Disaster Risk Reduction*, 141, 106172.

## See also

\[hm_empirical_correlogram()\], \[hm_simulate_sf_kle()\]

## Examples

``` r
if (FALSE) { # \dontrun{
corg  <- hm_empirical_correlogram(tau, D)
model <- hm_fit_correlation_model(corg)
model$best_family
model$all_fits
} # }
```
