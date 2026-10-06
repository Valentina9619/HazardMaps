# Plot frequency of best-fitting marginal distribution families

Draws a bar chart showing how many grid cells were assigned each
distribution family by \[hm_fit_marginals()\]. Reproduces the style of
Fig. A1 in Clavijo Mesa et al. (2026, Appendix A).

## Usage

``` r
hm_plot_dist_frequency(
  fits,
  title = NULL,
  bar_fill = NULL,
  bar_colour = "white",
  show_pct = TRUE
)
```

## Arguments

- fits:

  The data.frame returned by \[hm_fit_marginals()\].

- title:

  Optional plot title.

- bar_fill:

  Bar fill colour. If \`NULL\` (default), each bar is coloured by
  distribution family using the same palette as \[hm_plot_best_dist()\].

- bar_colour:

  Bar border colour.

- show_pct:

  Logical. If \`TRUE\`, adds percentage labels above each bar.

## Value

A \`ggplot\` object.

## See also

\[hm_fit_marginals()\], \[hm_plot_best_dist()\],
\[hm_plot_marginal_fit()\]

## Examples

``` r
if (FALSE) { # \dontrun{
fits <- hm_fit_marginals(X)
hm_plot_dist_frequency(fits)
} # }
```
