# Plot marginal fit diagnostics for a single grid cell

Produces a three-panel diagnostic figure for the best-fitting marginal
distribution at a chosen grid cell: (1) a histogram of the observed
values with the fitted PDF overlaid, (2) the empirical vs fitted CDF,
and (3) a quantile-quantile plot. Reproduces the style of Fig. 7 and
Appendix B (Fig. B1) in Clavijo Mesa et al. (2026).

## Usage

``` r
hm_plot_marginal_fit(
  X,
  fits,
  cell_id,
  title = NULL,
  x = NULL,
  n_bins = 30,
  colour_fit = "firebrick"
)
```

## Arguments

- X:

  A numeric matrix (same one passed to \[hm_fit_marginals()\]).

- fits:

  The data.frame returned by \[hm_fit_marginals()\].

- cell_id:

  Integer. The cell to diagnose (column index of `X`).

- title:

  Optional overall plot title. Defaults to the cell coordinates when
  \`x\` is supplied, or to the cell index otherwise.

- x:

  Optional \`hm_hazard\` object after \[hm_standardize_coords()\]. When
  supplied, the cell's geographic coordinates are shown in the title.

- n_bins:

  Number of histogram bins.

- colour_fit:

  Colour used for the fitted PDF and CDF lines.

## Value

A \`ggplot\` object (three panels arranged with
\[ggplot2::facet_wrap()\] on a long-format data frame).

## See also

\[hm_fit_marginals()\], \[hm_marginal_gof()\], \[hm_plot_best_dist()\],
\[hm_plot_dist_frequency()\]

## Examples

``` r
if (FALSE) { # \dontrun{
fits <- hm_fit_marginals(X)
hm_plot_marginal_fit(X, fits, cell_id = 1, x = x_ext)
} # }
```
