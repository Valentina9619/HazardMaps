# Compute spatial cluster sizes of extreme events

For each time step, counts the land cells where the hazard value meets
or exceeds `threshold`. This is the cluster-size metric \\B(t)\\ (Eq.
21, Appendix C of Clavijo Mesa et al., 2026) used to validate whether a
dependence model reproduces the observed spatial footprint.

## Usage

``` r
hm_cluster_size(X, threshold, time = NULL)
```

## Arguments

- X:

  Numeric matrix (\\T \times N\\) from \[hm_filter_land_cells()\].

- threshold:

  A single finite numeric value.

- time:

  Optional date vector of length \\T\\. Defaults to row indices.

## Value

A data.frame with columns `date`, `cluster_size`, and
`cluster_fraction`.

## References

Clavijo Mesa, M.V., Di Maio, F. & Zio, E. (2026). *Reliability
Engineering & System Safety*, 275, 112830.

## See also

\[hm_filter_land_cells()\], \[hm_select_extreme_events()\]

## Examples

``` r
if (FALSE) { # \dontrun{
X_land <- hm_filter_land_cells(X)
cs     <- hm_cluster_size(X_land, threshold = 25,
                           time = attr(X_land, "time"))
summary(cs$cluster_size)
} # }
```
