# Compute the empirical pairwise Kendall correlation matrix

Estimates Kendall's \\\hat{\tau}\\ between every pair of grid cells from
the Gaussian scores matrix `Z` from \[hm_pit_transform()\]. Implements
Eq. 20 of Clavijo Mesa et al. (2026, Copernicus paper).

## Usage

``` r
hm_empirical_kendall(Z, cell_idx = NULL)
```

## Arguments

- Z:

  Numeric matrix (\\T \times N\\) of centred Gaussian scores from
  \[hm_pit_transform()\].

- cell_idx:

  Optional integer vector of column indices. Defaults to all columns.

## Value

A symmetric \\N \times N\\ matrix of Kendall \\\hat{\tau}\\ values in
\\\[-1, 1\]\\.

## Details

For large grids (N \> 1000) supply a subsample via `cell_idx` to keep
computation feasible.

## References

Clavijo Mesa, M.V., Di Maio, F. & Zio, E. (2026). *Reliability
Engineering & System Safety*, 275, 112830.

## See also

\[hm_pit_transform()\], \[hm_empirical_correlogram()\]

## Examples

``` r
if (FALSE) { # \dontrun{
fits <- hm_fit_marginals(X)
pit  <- hm_pit_transform(X, fits)
tau  <- hm_empirical_kendall(pit$Z, cell_idx = 1:30)
range(tau)
} # }
```
