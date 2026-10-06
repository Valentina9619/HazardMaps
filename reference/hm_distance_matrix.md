# Compute pairwise haversine distances between grid cells

Returns a symmetric \\N \times N\\ distance matrix of great-circle
distances in kilometres between all pairs of land cells (Eq. 5 of
Clavijo Mesa et al., 2026, SF paper).

## Usage

``` r
hm_distance_matrix(coords)
```

## Arguments

- coords:

  A data.frame with `lat` and `lon` columns in decimal degrees.
  Typically `attr(X_land, "coords")` from \[hm_filter_land_cells()\].

## Value

A symmetric numeric matrix with diagonal zero.

## References

Clavijo Mesa, M.V., Broggi, M., Di Maio, F. & Zio, E. (2026).
*International Journal of Disaster Risk Reduction*, 141, 106172.

## See also

\[hm_empirical_correlogram()\], \[hm_filter_land_cells()\]

## Examples

``` r
if (FALSE) { # \dontrun{
X_land <- hm_filter_land_cells(X)
D      <- hm_distance_matrix(attr(X_land, "coords"))
range(D)
} # }
```
