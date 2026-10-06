# Plot the empirical spatial correlogram

Draws the empirical correlogram produced by
\[hm_empirical_correlogram()\]: average Kendall \\\hat{\tau}\\ as a
function of inter-site distance, with point size proportional to the
number of cell pairs in each bin. Reproduces the style of Fig. 8 in
Clavijo Mesa et al. (2026, SF paper).

## Usage

``` r
hm_plot_correlogram(
  corg,
  title = NULL,
  colour = "#2166AC",
  show_points = TRUE,
  show_zero = TRUE,
  ylim = NULL
)
```

## Arguments

- corg:

  A data.frame returned by \[hm_empirical_correlogram()\], with columns
  \`bin_center\`, \`mean_tau\`, and \`n_pairs\`.

- title:

  Optional plot title.

- colour:

  Line and point colour. Defaults to \`"#2166AC"\` (steel blue).

- show_points:

  Logical. If \`TRUE\` (default), draws points at each bin centre sized
  by the number of pairs.

- show_zero:

  Logical. If \`TRUE\` (default), adds a horizontal dashed line at
  \\\tau = 0\\.

- ylim:

  Optional numeric vector of length 2 for the y-axis range. If \`NULL\`
  (default), the range is set automatically.

## Value

A \`ggplot\` object.

## Details

A horizontal dashed reference line is drawn at \\\tau = 0\\ to help
identify the distance at which spatial dependence becomes negligible.
Point size is scaled by \`sqrt(n_pairs)\` so bins with more pairs are
visually more prominent without dominating the plot.

## References

Clavijo Mesa, M.V., Broggi, M., Di Maio, F. & Zio, E. (2026).
Inoperability assessment of interdependent critical infrastructures
exposed to natural hazards considering climate change. *International
Journal of Disaster Risk Reduction*, 141, 106172.

## See also

\[hm_empirical_correlogram()\], \[hm_empirical_kendall()\],
\[hm_distance_matrix()\]

## Examples

``` r
if (FALSE) { # \dontrun{
X_land <- hm_filter_land_cells(X)
fits   <- hm_fit_marginals(X_land)
pit    <- hm_pit_transform(X_land, fits)
tau    <- hm_empirical_kendall(pit$Z)
D      <- hm_distance_matrix(attr(X_land, "coords"))
corg   <- hm_empirical_correlogram(tau, D)

hm_plot_correlogram(corg,
  title = "Empirical correlogram — Copernicus wind gust")
} # }
```
