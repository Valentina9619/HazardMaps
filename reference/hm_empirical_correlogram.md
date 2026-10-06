# Compute the empirical spatial correlogram

Bins pairwise Kendall correlations by inter-site distance to give
average correlation as a function of distance. The result is the key
input for fitting a covariance model in \[hm_fit_correlation_model()\],
following Eq. 5–6 of Clavijo Mesa et al. (2026, SF paper).

## Usage

``` r
hm_empirical_correlogram(tau, D, n_bins = 30L, max_dist = NULL)
```

## Arguments

- tau:

  Symmetric numeric matrix of pairwise Kendall correlations from
  \[hm_empirical_kendall()\].

- D:

  Symmetric numeric distance matrix in km from \[hm_distance_matrix()\].
  Must have the same dimensions as `tau`.

- n_bins:

  Number of equal-width distance bins. Defaults to 30.

- max_dist:

  Maximum distance in km to include. Defaults to the observed maximum.

## Value

A data.frame with columns `bin_center`, `mean_tau`, `n_pairs`, `bin_lo`,
and `bin_hi`.

## References

Clavijo Mesa, M.V., Broggi, M., Di Maio, F. & Zio, E. (2026).
*International Journal of Disaster Risk Reduction*, 141, 106172.

## See also

\[hm_empirical_kendall()\], \[hm_distance_matrix()\],
\[hm_fit_correlation_model()\]

## Examples

``` r
if (FALSE) { # \dontrun{
tau  <- hm_empirical_kendall(pit$Z, cell_idx = 1:50)
D    <- hm_distance_matrix(coords[1:50, ])
corg <- hm_empirical_correlogram(tau, D, n_bins = 20)
head(corg)
} # }
```
