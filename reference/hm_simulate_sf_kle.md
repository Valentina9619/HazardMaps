# Simulate hazard scenarios via Karhunen-Loève Expansion

Generates `n_scenarios` synthetic hazard fields from the KLE of the
fitted covariance model, following Eq. 8–11 of Clavijo Mesa et al.
(2026, SF paper): build \\\bar{C}\\, eigen-decompose and truncate to
`var_retain` variance, sample KLE coefficients, back-transform via
\\\Phi\\ and the inverse marginal CDFs.

## Usage

``` r
hm_simulate_sf_kle(
  corrmodel,
  D,
  fits,
  n_scenarios = 1000L,
  var_retain = 0.95,
  seed = NULL
)
```

## Arguments

- corrmodel:

  An `hm_corrmodel` from \[hm_fit_correlation_model()\].

- D:

  Symmetric distance matrix in km from \[hm_distance_matrix()\]. Must
  have `nrow(fits)` rows.

- fits:

  data.frame from \[hm_fit_marginals()\] for the same land cells.

- n_scenarios:

  Number of synthetic scenarios. Defaults to 1000.

- var_retain:

  Variance fraction retained in KLE truncation. Defaults to 0.95.

- seed:

  Optional integer random seed.

## Value

A numeric matrix (\\S \times N\\) of simulated hazard intensities.
Attributes `family`, `k_modes`, `var_retain`, and `var_explained` are
attached.

## References

Clavijo Mesa, M.V., Broggi, M., Di Maio, F. & Zio, E. (2026).
*International Journal of Disaster Risk Reduction*, 141, 106172.

## See also

\[hm_fit_correlation_model()\], \[hm_fit_marginals()\]

## Examples

``` r
if (FALSE) { # \dontrun{
model     <- hm_fit_correlation_model(corg)
scenarios <- hm_simulate_sf_kle(model, D, fits,
                                 n_scenarios = 1000, seed = 42)
attr(scenarios, "k_modes")
attr(scenarios, "var_explained")
} # }
```
