# Compute the empirical spatial correlogram

Bins pairwise Kendall correlations by inter-site distance to give
average correlation as a function of distance. The result is the key
input for fitting a covariance model in \[hm_fit_correlation_model()\].

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
f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
                 package = "HazardMaps")
x      <- hm_read_netcdf(f, var = "i10fg")
x      <- hm_standardize_coords(x)
x      <- hm_decode_time(x)
x_ext  <- hm_select_extreme_events(x, threshold = 25)
X      <- hm_to_matrix(x_ext)
X_land <- hm_filter_land_cells(X)
fits   <- hm_fit_marginals(X_land)
pit    <- hm_pit_transform(X_land, fits)
tau    <- hm_empirical_kendall(pit$Z)
D      <- hm_distance_matrix(attr(X_land, "coords"))
corg   <- hm_empirical_correlogram(tau, D, n_bins = 20L)
head(corg)
} # }
```
